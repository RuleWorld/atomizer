import readBNGXML
import progressbar
import seaborn as sns
import matplotlib.pyplot as plt
import os
import fnmatch
import pandas
from scipy.stats import kendalltau
from scipy import stats
import argparse

xkcd = False

def getModelStructures(bngxml):
    structures = readBNGXML.parseFullXML(bngxml)
    return structures

def constructHistogram(data, fileName, xlabel, ylabel, bins=10):
    """
    constructs a histogram based on the information in data
    """
    _, axs = plt.subplots(1, 1, sharex=True, figsize=(8, 6))

    plt.clf()
    #axs.set_xlim(xmin=0,xmax=bins[-1])
    sns.distplot(data, kde=False, rug=True, bins=bins, color="g")
    #plt.hist(ratomization)
    
    plt.savefig(fileName)

def categorizeStatistics(effort):
    if effort == 0:
        return '0'
    elif effort >=1 and effort <= 5:
        return '1-5'
    elif effort >5 and effort <= 10:
        return '6-10'
    elif effort >10:
        return '>10'

def create2DdensityPlot(dataframe, columns, axisNames, outputfilename, plotType=sns.kdeplot, xlim=(0, 1), ylim=(0, 1)):
    """
    creates a 2d density plot given a dataframe and two columns.
    saves the image in <outputfilename>
    """
    plt.clf()
    _, _ = plt.subplots(1, 1, sharex=True, figsize=(8, 6))
    g = sns.JointGrid(columns[0], columns[1], dataframe, space=0)
    g.plot_marginals(sns.distplot, color="g", bins=None)
    g.plot_joint(plotType, cmap="Greens", shade=True, n_levels=20)

    g.set_axis_labels(axisNames[0], axisNames[1])

    #ax = g.ax_joint
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #g.ax_marg_x.set_xscale('log')
    #g.ax_marg_y.set_yscale('log')
    g.annotate(stats.pearsonr)
    plt.savefig(outputfilename)



def getActiveMolecules(rules):
    activeMolecules = set([])
    for rule in rules:
        reactants = set([])
        products = set([])
        reactants = set([x.name for y in rule[0].reactants for x in y.molecules if x.name != '0'])
        products = set([x.name for y in rule[0].reactants for x in y.molecules if x.name != '0'])
        # stuff that just moves around doesn't count
        actions =  [x.action for x in rule[0].actions if x.action != 'ChangeCompartment']
        if reactants and products and len(actions) > 0:
            activeMolecules = activeMolecules | reactants
            activeMolecules = activeMolecules | products

    return activeMolecules

def underscoreTest(moleculeTypes, activeMolecules):
    '''
    notice how the method's name doesn't have an underscore.
    '''
    totalMolecules = len(activeMolecules)
    underscoredMolecules = 0
    for moleculeType in moleculeTypes:
        if moleculeType.name not in activeMolecules:
            continue
        tokens = moleculeType.name.split('_')
        tokens = [x for x in tokens if len(x) > 0]
        if len(tokens) > 3:
            underscoredMolecules += 1
    return underscoredMolecules, totalMolecules

def getFiles(directory, extension):
    """
    Gets a list of <*.extension> files. include subdirectories and return the absolute 
    path. also sorts by size.
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.{0}'.format(extension)):
            matches.append([os.path.join(os.path.abspath(root), filename),os.path.getsize(os.path.join(root, filename))])

    #sort by size
    matches.sort(key=lambda filename: filename[1], reverse=False)
    
    matches = [x[0] for x in matches]

    return matches



def batchComparison(directory):
    normalizedScore = {}
    score = {}
    categorizedScore = {}
    xmlfiles = getFiles(directory, 'xml')
    progress = progressbar.ProgressBar()

    atomizationDB = pandas.read_hdf('{0}DB.h5'.format(directory), 'atomization')

    for xmlidx in progress(range(len(xmlfiles))):

        xml = xmlfiles[xmlidx]
        filename = os.sep.join(xml.split(os.sep)[-2:])

        if  filename not in atomizationDB.index:
            continue
        xmlstructures = getModelStructures(xml)

        if  all([len(x[0].split('_')) <= 2 for x in xmlstructures['observables']]):
            continue

        activeMolecules = getActiveMolecules(xmlstructures['rules'])
        if len(activeMolecules) > 0:
            modelscore, total = underscoreTest(xmlstructures['molecules'], activeMolecules)
            score[filename] = modelscore
            normalizedScore[filename] = modelscore*1.0/total
            categorizedScore[filename] = categorizeStatistics(modelscore)

    atomizationDB['nonatoscore'] = pandas.Series(score, index=atomizationDB.index)
    atomizationDB['categorized_nonatoscore'] = pandas.Series(categorizedScore, index=atomizationDB.index)
    atomizationDB['normalized_nonatoscore'] = pandas.Series(normalizedScore, index=atomizationDB.index)
    atomizationDB.to_hdf('{0}DB.h5'.format(directory),'atomization')
    return atomizationDB

def createPlots(atomizationDB):

    plt.clf()

    #sns.set(font_scale=1.5) 
    sns.set_context("paper", font_scale=3)
    sns.set_style("ticks")
    if xkcd:
        plt.xkcd()

    prunnedDB =  atomizationDB[atomizationDB.nonatoscore.notnull()]

    g = sns.factorplot(x="categorized_nonatoscore", data=prunnedDB, kind="count",
                   palette="BuGn_d", size=7, aspect=0.9, order=['0', '1-5' , '6-10', '>10'])

    plt.xlabel("Unatomized species\n ({0} models)".format(len(prunnedDB)), fontsize=32,fontweight='bold')
    plt.ylabel("Number of Models", fontsize=32,fontweight='bold')
    #g.set_axis_labels("Number of unatomized species\n in models with complex formation ({0} models)".format(len(prunnedDB)), "Number of Models")

    sns.despine()
    if xkcd:
        g.fig.savefig('underscorebarplot_xkcd.png',bbox_inches='tight')
    else:
        g.fig.savefig('underscorebarplot.png',bbox_inches='tight')
    #nonatoscore = [x for x,y in zip(atomizationDB['nonatoscore'], atomizationDB['nonatoscore'].notnull()) if y]
    #normalized_nonatoscore = [x for x,y in zip(atomizationDB['normalized_nonatoscore'], atomizationDB['normalized_nonatoscore'].notnull()) if y]

    #constructHistogram(nonatoscore, 'unatomizedScore.png', 'Non-atomized ratio','Number of models')
    #constructHistogram(normalized_nonatoscore, 'unatomizedMolecules.png', 'Non-atomized molecules','Number of models')

    #simple scatterplot comparing underscore to effort

    #plt.clf()
    #_, axs = plt.subplots(1, 1, sharex=True, figsize=(8, 6))

    #shifted2DDensityPlot(prunnedDB['nonatoscore'], prunnedDB['effort'], [5,10,15,20], [5,10,15,20], 'nonatoscorevseffort.png')
    #create2DdensityPlot(prunnedDB, ['nonatoscore', 'effort'], ['number of unatomized species in models with complex formation (that use the underscore convention) score', 'Amount of user information required by atomizer'], 'nonatoscorevseffort.png')


    #sns.lmplot("nonatoscore", "y", data=atomizationDB, hue='categorized_effort', fit_reg=False)
def defineConsole():
    parser = argparse.ArgumentParser(description='SBML to BNGL translator')
    parser.add_argument('-i','--input',type=str,help='input SBML model1',required=True)
    parser.add_argument('-x','--xkcd',action='store_true',help='xkcd type plots')
    #parser.add_argument('-o','--output-file',type=str,help='output SBML file',required=True)
    return parser    


if __name__ == "__main__":
    #xmlstructures = getModelStructures('curated/BIOMD0000000470.xml.xml')
    #activeMolecules = getActiveMolecules(xmlstructures['rules'])
    parser = defineConsole()
    namespace = parser.parse_args()
    xkcd = namespace.xkcd
    directory = 'curated'
    #atomizationDB = batchComparison(directory)
    atomizationDB = pandas.read_hdf('{0}DB.h5'.format(directory), 'atomization')

    createPlots(atomizationDB)
    #print atomizationDB.query('(nonatoscore > 50)')
