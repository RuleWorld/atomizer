import pandas
from scipy import stats
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
import progressbar
from richContactMap import reactionBasedAtomization, stoichiometryAnalysis, extractActiveMolecules, getValidFiles
import readBNGXML
from collections import Counter
import SBMLparser.utils.annotationExtractor as annEx
from scipy.stats import kendalltau

sns.set_style("white")

def constructHistogram(data, fileName, xlabel, ylabel, bins=10):
    """
    constructs a histogram based on the information in data
    """
    _, axs = plt.subplots(1, 1, sharex=True, figsize=(8, 6))

    plt.clf()
    sns.set_palette("BuGn_d")
    if type(bins) != int:
        axs.set_xlim(xmin=0,xmax=bins[-1])
    sns.distplot(data, kde=False, rug=False, bins=bins, hist_kws=dict(alpha=1))
    #plt.hist(ratomization)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    
    plt.savefig(fileName)


def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)

def create2DdensityPlot(dataframe, columns, axisNames, outputfilename, plotType=sns.kdeplot, xlim=(0, 1), ylim=(0, 1)):
    """
    creates a 2d density plot given a dataframe and two columns.
    saves the image in <outputfilename>
    """
    plt.clf()
    _, _ = plt.subplots(1, 1, sharex=True, figsize=(8, 6))
    g = sns.JointGrid(columns[0], columns[1], dataframe, xlim=xlim, ylim=ylim, space=0)
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


def createHexBin(dataframe,columns,axisnames,outputfilename,xlim=(0,1),ylim=(0,1)):
    plt.clf()

    g = sns.JointGrid(columns[0], columns[1], dataframe, space=0)
    g.ax_marg_x.hist(dataframe[columns[0]], bins=np.arange(xlim[0], xlim[1]))
    g.ax_marg_y.hist(dataframe[columns[1]], bins=np.arange(ylim[0], ylim[1], orientation="horizontal"))

    #g.ax_marg_x.hist(x, bins=np.arange(0, 60)
    #g.ax_marg_y.hist(y, bins=np.arange(0, 1000, 10), orientation="horizontal")
    g.plot_joint(plt.hexbin, gridsize=25, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], cmap="Blues")
    #g.fig.savefig("/Users/mwaskom/Desktop/jointgrid.png", bbox_inches="tight")

    
    #f, _ = plt.subplots(1, 1, sharex=True, figsize=(8, 6))
    #g = sns.JointGrid(columns[0], columns[1], dataframe, xlim=xlim, ylim=ylim, space = 0)
    #g.plot_marginals(sns.distplot, color="g")
    #g.plot_joint(plt.hexbin, cmap="Greens", extent=[0, np.max(dataframe[columns[0]]), 0, np.max(dataframe[columns[1]])])
    #g.annotate(stats.pearsonr)
    #sns.jointplot(dataframe[columns[0]], dataframe[columns[1]], kind="hex", stat_func=stats.pearsonr, color="g", gridsize=8)
    g.fig.savefig(outputfilename)

def create1Ddensityplot(data, outputfilename):
    plt.clf()
    f, (ax1) = plt.subplots(1, 1, sharex=True, figsize=(8, 6))
    # with sns.axes_style("white"):
    #sns.jointplot("compression", "wiener index",atomizationInfo, kind="kde");
    sns.kdeplot(data, shade=True, ax=ax1, clip=(0, 1), bw=0.5)
    plt.savefig(outputfilename)




def reactionBasedAtomizationDistro(directory):
    '''
    calculates a rection atomization based metric:
    ration of atomized reactions (non syndeg) in a model
    '''
    syndelArray = []
    atomizedDistro = []
    nonAtomizedDistro = []
    atomizationDB = pandas.DataFrame()


    # generate bng-xml
    # generateBNGXML(directory)

    print 'reading bng-xml files'
    xmlFiles = getValidFiles(directory, 'xml')

    print 'analyzing {0} bng-xml files'.format(len(xmlFiles))
    progress = progressbar.ProgressBar()
    validFiles = 0
    for i in progress(range(len(xmlFiles))):

        xml = xmlFiles[i]
    # for xml in xmlFiles:
        try:
            # console.bngl2xml('complex/output{0}.bngl'.format(element),timeout=10)
            try:

                structures = readBNGXML.parseFullXML(xml)
                rules = structures['rules']
                observables = structures['observables']
                molecules = structures['molecules']
            except IOError:
                print xml
                continue
            atomizedProcesses, weight = reactionBasedAtomization(rules)
            ato, nonato = stoichiometryAnalysis(rules)
            atomizedDistro.extend(ato)
            nonAtomizedDistro.extend(nonato)
            # if (2,1) in nonato:
            #    interesting.append(element)
            score = atomizedProcesses * 1.0 / weight if weight != 0 else 0
            #totalRatomizedProcesses += atomizedProcesses
            #totalReactions += len(rules)
            #totalProcesses += weight

            # calculate yield
            activeMolecules = extractActiveMolecules(rules)
            activeMoleculeTypes = [x for x in molecules if x.name in activeMolecules]
            yieldValue = len([x for x in activeMoleculeTypes if len(
                x.components) > 0]) * 1.0 / len(activeMoleculeTypes) if len(activeMoleculeTypes) > 0 else 0

            # syndel value
            syndelValue = 1 - (len(rules) - weight) * 1.0 / len(rules) if len(rules) > 0 else 0

            atomizationDB.set_value(xml, 'score', score)
            atomizationDB.set_value(xml, 'weight', weight)
            atomizationDB.set_value(xml, 'length', len(rules))
            atomizationDB.set_value(xml, 'yild', yieldValue)
            atomizationDB.set_value(xml, 'syndel', syndelValue)
            atomizationDB.set_value(xml, 'numspecies', len(observables))
            validFiles += 1
        except IOError:
            print 'io'
            continue
    print 'found {0} models i could extract info from'.format(validFiles)

    return atomizationDB


def extractAnnotationsFromModelSet(modelList):
    modelAnnotationsCounter = Counter()
    for model in modelList:
        annotationExtractor = annEx.AnnotationExtractor(model)
        modelAnnotations = annotationExtractor.getModelAnnotations()
        
        speciesAnnotations = annotationExtractor.getAnnotationSystem()
        #print speciesAnnotations
        speciesAnnotations = set([z for x in speciesAnnotations for y in speciesAnnotations[x] for z in speciesAnnotations[x][y]])
        modelAnnotationsCounter.update(speciesAnnotations)

    print modelAnnotationsCounter.most_common(20)



def constructPlots(atomizationDB):
    """
    Given a pandas data frame object it creates a series of histogram and kde plots describing the characteristics of a set of atomized
    models
    """
    constructHistogram(atomizationDB['syndel'], '{0}/syndelHist.png'.format(outputDir), 'Fraction of non syn-del reactions', 'Number of models')             
    constructHistogram(atomizationDB['yild'], '{0}/yieldHist.png'.format(outputDir), 'Yield score', 'Number of models')             
    constructHistogram(atomizationDB['score'], '{0}/atomizationHist.png'.format(outputDir), 'Percentage of reactions with mechanistic processes', 'Number of models')             
    create2DdensityPlot(atomizationDB, ['score', 'yild'], ['Atomization score', 'Yield score'], '{0}/atomizationvsyield.png'.format(outputDir))
    create2DdensityPlot(atomizationDB, ['syndel', 'yild'], ['Percentage or non-syndel reactions', 'Yield score'], '{0}/syndelvsyield.png'.format(outputDir))

    #createHexBin(atomizationDB, ['syndel', 'yild'], ['Percentage or non-syndel reactions', 'Yield score'], '{0}/syndelvsyieldhex.png'.format(outputDir))

if __name__ == "__main__":
    folder = 'curated'
    # calculate atomization information
    atomizationDB = reactionBasedAtomizationDistro(folder)
    atomizationDB.to_hdf('{0}DB.h5'.format(folder),'atomization')

    outputDir = 'testCurated'
    # read info
    #atomizationDB = pandas.read_hdf('{0}DB.h5'.format(folder), 'atomization')

    
    # construct plots
    #constructPlots(atomizationDB)

    #testSet = list(atomizationDB.query('(yild < 0.6) & (syndel > 0.8)').index)
    #testSet = ['XMLExamples/{0}'.format(x[:-4]) for x in testSet]
    #print extractAnnotationsFromModelSet(testSet)

