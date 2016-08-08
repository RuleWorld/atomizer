import os.path
import collections
import cPickle as pickle
from copy import copy
import math
import progressbar
import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import spectral_clustering
import difflib
import fnmatch
import matplotlib.pyplot as plt
import pandas
import seaborn as sns
sns.set_style("ticks")
sns.despine()
sns.set_context("paper", font_scale=2, rc={"lines.linewidth": 2.5})
sns.set_style("ticks")
sns.despine()
import sys
import os
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), 'SBMLparser'))
from analyzeModelSet import reactionBasedAtomizationFile as ratofile
import SBMLparser.utils.readBNGXML as readBNGXML
import SBMLparser.utils.smallStructures as smallStructures

def extractMoleculeTypesFromFile(fileName):
    species, rules, par = readBNGXML.parseXML(fileName)
    return species

def getValidFiles(directory, extension):
    """
    Gets a list of bngl files that could be correctly translated in a given 'directory'
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.{0}'.format(extension)):
            matches.append(os.path.join(root, filename))
    return matches

def extractRulesFromFile(fileName):
    species, rules, par = readBNGXML.parseXML(fileName)
    return rules

def matchAnnotationsToSpecies(moleculeArray, annotationDictionary):
    '''
    Given a list of molecule types and a set of annotations, this method
    will match and create a dictionary containing  the annotations for
    the intersection
    '''
    parsedDictionary = {}
    finalDictionary = {}
    parsedDictionary = copy(annotationDictionary)
    for annotationName in annotationDictionary:
        #speciesName = '_'.join(annotationName.split('_')[0:-1])
        #speciesId = annotationName.split('_')[-1]
        parsedDictionary[speciesName] = annotationDictionary[annotationName]
        parsedDictionary[speciesId] = annotationDictionary[annotationName]

    # print parsedDictionary
    for molecule in moleculeArray:
        if molecule.name in annotationDictionary:
            finalDictionary[molecule.name] = parsedDictionary[molecule.name]

    return parsedDictionary


def createNeighborhoodDictionary(moleculeArray):
    neighborhoodDictionary = collections.defaultdict(list)
    moleculeNames = {x.name.lower(): x.name for x in moleculeArray}

    for molecule in moleculeArray:
        for component in molecule.components:
            if component.name in moleculeNames:
                neighborhoodDictionary[molecule.name].append(
                    moleculeNames[component.name])
    return neighborhoodDictionary


def resolveAnnotation(annotations):
    '''
    
    '''
    with open('parsedAnnotations.dump', 'rb') as f:
        parsedAnnotations = pickle.load(f)

    newAnnotations = {}
    for annotationKey in annotations:
        for annotation in annotations[annotationKey]:
            if annotation in parsedAnnotations and \
                    parsedAnnotations[annotation] != []:
                newAnnotations[annotationKey] = parsedAnnotations[annotation]
    return newAnnotations


def extractMoleculeTypes(folderName, annotationsFolder, includeAnnotations=True):
    '''
    return a list of molecule types structures
    from a group of BNG-XML files in a folder <folderName>
    '''
    moleculeTypesArray = []
    progress = progressbar.ProgressBar()
    #with open('annotations.dump', 'rb') as f:
    #    annotationsArray = pickle.load(f)

    with open(os.path.join(annotationsFolder, 'annotationDictionary.dump'), 'rb') as f:
        annotationsArray = pickle.load(f)
    bngxmlFiles = getValidFiles(folderName, 'xml')
    for element in progress(range(0, len(bngxmlFiles))):
        fileName = bngxmlFiles[element]
        try:
        
            moleculeTypes = extractMoleculeTypesFromFile(fileName)
            observablesLen = readBNGXML.getNumObservablesXML(fileName)
            ratoscore = ratofile(fileName, None, None)[1]
            if includeAnnotations:
                #annotations = matchAnnotationsToSpecies(moleculeTypes, 
                #                                        annotationsArray[fileName])
                
                try:
                    annotations = annotationsArray[os.path.join(annotationsFolder, fileName[:-4].split('/')[1])]
                    resolvedAnnotations = resolveAnnotation(annotations)
                except KeyError:
                    print os.path.join([annotationsFolder, fileName[:-4].split('/')[1]])
                    continue
                    

            else:
                annotations = []
                resolvedAnnotations = []
            neighborhoodDictionary = createNeighborhoodDictionary(
                moleculeTypes)
            moleculeTypesArray.append(
                [moleculeTypes, annotations, resolvedAnnotations, 
                    neighborhoodDictionary, fileName, observablesLen, ratoscore])
        except IOError:
            continue
    return moleculeTypesArray

from pandas import Series

def extractProcessDistribution(folderName, cluster=False):
    '''
    return a list of molecule types structures
    from a group of BNG-XML files in a folder <folderName>
    '''
    processArray = []
    progress = progressbar.ProgressBar()
    #with open('annotations.dump', 'rb') as f:
    #    annotationsArray = pickle.load(f)

    #with open(os.path.join(annotationsFolder, 'annotationDictionary.dump'), 'rb') as f:
    #    annotationsArray = pickle.load(f)

    bngxmlFiles = getValidFiles(folderName[0], 'xml')
    modelHistogram = collections.defaultdict(list)
    actionDict = collections.defaultdict(list)
    for element in progress(range(0, len(bngxmlFiles))):
        if not cluster:
            actionHistogram = {'DeleteBond': 0, 'Add': 0, 'StateChange': 0, 
            'AddBond': 0, 'Delete': 0, 'ChangeCompartment':0}
        else:
            actionHistogram = {'Add/Delete': 0, '(Add/Delete)Bond':0, 'StateChange':0, 'ChangeCompartment':0}
        fileName = bngxmlFiles[element]
        try:
            rules = extractRulesFromFile(fileName)
            for rule in rules:
                for action in rule[0].actions:
                    if not cluster:
                        actionHistogram[action.action] += 1
                    else:
                        if action.action in ['Add', 'Delete']:
                            actionHistogram['Add/Delete'] += 1
                        elif action.action in ['AddBond', 'DeleteBond']:
                            actionHistogram['(Add/Delete)Bond'] += 1
                        else:
                            actionHistogram[action.action] += 1
            del actionHistogram['ChangeCompartment']
            totalProcesses = sum([actionHistogram[x] for x in actionHistogram])
            if totalProcesses == 0:
                continue
            actionHistogram = {x:actionHistogram[x]*1.0/totalProcesses for x in actionHistogram}
            for x in actionHistogram:
                actionDict['process'].append(x)
                actionDict['fraction'].append(actionHistogram[x])
                actionDict['database'].append(folderName[1])
                actionDict['file'].append(fileName)
            #actionHistogram.sort(key=lambda x:x[0])
            for element in actionHistogram:
                modelHistogram[element].append(actionHistogram[element])
            
        except IOError:
            continue


    
    for x in modelHistogram:
        modelHistogram[x] = [np.mean(modelHistogram[x]), np.std(modelHistogram[x])]

    return modelHistogram, actionDict
    '''
    s = Series([x for x in actionHistogram.elements()])
    vc = s.value_counts()
    ax = vc.plot(kind='bar', fontsize=18, rot=0)
    fig = ax.get_figure()
    fig.set_size_inches(18.5, 12.0)
    fig.savefig('{0}process.png'.format(folderName))
    print '{0}process.png'.format(folderName)
    #print actionHistogram
    '''

def compareStrings(name1, name2, minThreshold=0.8):
    if len(name1) > 3 and len(name2) > 3:
        similarity = difflib.SequenceMatcher(a=name1, b=name2).quick_ratio()
        #if similarity < minThreshold:
        #    similarity = -1
        affinity =  1000 * similarity
        return affinity

    return -1e3


def compareAnnotations(annotation1, annotation2):
    affinity = 0
    for ann1 in annotation1:
        for ann2 in annotation2:
            affinity += compareStrings(ann1, ann2)
    return affinity


def distanceFunction(moleculeA, moleculeB):
    moleculeAStructure = moleculeA[0]
    moleculeBStructure = moleculeB[0]
    moleculeAAnn = moleculeA[1]
    moleculeBAnn = moleculeB[1]
    moleculeANeigh = moleculeA[3]
    moleculeBNeigh = moleculeB[3]


    affinity = compareStrings(moleculeAStructure, moleculeBStructure)
    
    affinity += 0.05 * compareAnnotations(moleculeAAnn, moleculeBAnn)
    for neighborA in moleculeANeigh:
        for neighborB in moleculeBNeigh:
            affinity += 0.1 * compareStrings(neighborA, neighborB, 0.7)

    return affinity


def createSimilarityMatrix(modelMoleculeTypeList):
    similarityMatrix = np.zeros(
        (len(modelMoleculeTypeList), len(modelMoleculeTypeList)))

    progress = progressbar.ProgressBar()
    for index1 in progress(range(0, len(modelMoleculeTypeList))):
        similarityMatrix[index1][index1] = 1000
        for index2 in range(index1+1, len(modelMoleculeTypeList)):
            similarityMatrix[index1][index2] = distanceFunction(modelMoleculeTypeList[index1], 
                                                                modelMoleculeTypeList[index2])
            similarityMatrix[index2][index1] = similarityMatrix[index1][index2]


    return similarityMatrix


def moleculeTypeNamesConsolidation(moleculetypestotalarray):
    modelMoleculeTypeList = []
    progress = progressbar.ProgressBar()

    for index in progress(range(0, len(moleculetypestotalarray))):
        model = moleculetypestotalarray[index]
        moleculeTypes = model[0]
        annotations = model[1]
        resolvedAnnotations = model[2]
        neighborhood = model[3]
        modelNumber = model[4]
        observablesLen = model[5]
        for molecule in moleculeTypes:
            if molecule.name in resolvedAnnotations:
                rann = resolvedAnnotations[molecule.name]
            else:
                rann = []
            if molecule.name in annotations:
                ann = annotations[molecule.name]
            else:
                ann = []
            modelMoleculeTypeList.append([molecule.name, 
                                          rann, modelNumber, 
                                          neighborhood[molecule.name], ann, observablesLen])
    return modelMoleculeTypeList


def rawAnnotationCoverage(directory):
    
    with open(os.path.join(directory, 'annotationDictionary.dump'), 'rb') as f:
        annotationInformation = pickle.load(f)
    print(sum([1 for x in annotationInformation.keys() for y in annotationInformation[x] if annotationInformation[x][y] != []]), sum([1
         for x in annotationInformation.keys() for y in annotationInformation[x]]))
    annotationCounter = collections.Counter()
    modelsWithoutAnnotation = collections.Counter()
    for model in annotationInformation:
        for species in annotationInformation[model]:
            if annotationInformation[model][species] == []:
                modelsWithoutAnnotation.update([model.split('/')[2]])
            for annotation in annotationInformation[model][species]:
                try:
                    annotationCounter.update([annotation.split('/')[3]])
                except:
                    continue
    print(annotationCounter)
    #print('no annotation', modelsWithoutAnnotation)


def preliminaryAnalysis(directory='new_non_curated', directory2=None):
    '''
    print('building data structures...')
    moleculeTypesArray = extractMoleculeTypes(directory, directory2, True)
    with open(os.path.join(directory, 'moleculeTypeDataSet.dump'), 'wb') as f:
        pickle.dump(moleculeTypesArray, f)
    '''
    with open(os.path.join(directory, 'moleculeTypeDataSet.dump'), 'rb') as f:
        moleculeTypesArray = pickle.load(f)

    #clusterMoleculeTypes(moleculeTypesArray)
    modelMoleculeTypeList = moleculeTypeNamesConsolidation(moleculeTypesArray)
    '''
    print('building similarity matrix...')
    similarityMatrix = createSimilarityMatrix(modelMoleculeTypeList)
    with open(os.path.join(directory, 'modelMoleculeTypeSimiliarityMatrix'), 'wb') as f:
        pickle.dump(similarityMatrix, f)
    
    with open(os.path.join(directory, 'modelMoleculeTypeSimiliarityMatrix'), 'rb') as f:
        similarityMatrix = pickle.load(f)
    
    print('running clustering algorithm...')
    af = AffinityPropagation(affinity="precomputed", preference=500).fit(similarityMatrix)
    #af = spectral_clustering(affinity="precomputed").fit(similarityMatrix)


    with open(os.path.join(directory, 'clusteringResult.dump'), 'wb') as f:
        pickle.dump(af, f)
    

    with open(os.path.join(directory, 'clusteringResult.dump'), 'rb') as f:
        af = pickle.load(f)
    cluster_centers_indices = af.cluster_centers_indices_

        
    print('final results:')
    labels = af.labels_

    n_clusters_ = len(cluster_centers_indices)

    print('Estimated number of clusters: %d' % n_clusters_)
    '''
    print('Original number of molecule types: {0}'.format(len(modelMoleculeTypeList)))
    print('Original number of species: {0}'.format(sum([x[5] for x in moleculeTypesArray])))
    print('Number of molecule types with annotations: {0}'.format(len([x for x in modelMoleculeTypeList if x[4] != []])))
    print('number of models analyzed: {0}'.format(len(moleculeTypesArray)))
    moleculeDictionary = collections.defaultdict(list)
    '''
    for element, label in zip(modelMoleculeTypeList, labels):
        moleculeDictionary[label].append(element)
    '''
    #print('molecule types without annotations: {0}'.format(([x for x in modelMoleculeTypeList if x[4] == []])))
    with open(os.path.join(directory, 'finalDictionary'), 'wb') as f:
        pickle.dump(moleculeDictionary, f)




def componentAnalysis(directory, atomizeThreshold=0):
    componentCount = []
    bindingCount = []
    stateCount = []
    atoarray = []
    with open(os.path.join(directory, 'moleculeTypeDataSet.dump'), 'rb') as f:
        moleculeTypesArray = pickle.load(f)
    for model in moleculeTypesArray:
        if model[-1] < atomizeThreshold:
            continue
        modelComponentCount = [len(x.components) for x in model[0]]
        bindingComponentCount = [len([y for y in x.components if len(y.states) == 0])
                                 for x in model[0]]

        modificationComponentCount = [sum([max(1, len(y.states)) for y in x.components])
                                      for x in model[0]]

        bindingCount.append(bindingComponentCount)
        stateCount.append(modificationComponentCount)
        componentCount.append(modelComponentCount)
        atoarray.append(model[-1])
    #print [(np.mean(x), np.std(x)) for x in componentCount]

    
    
    componentCount = np.array(componentCount)

    bindingCount = np.array(bindingCount)
    stateCount = np.array(stateCount)
    atoarray = np.array(atoarray)
    #interestingCount = [index for index, x in enumerate(componentCount) if np.mean(x) >= atomizeThreshold]
    #componentCount = componentCount[interestingCount]
    #bindingCount = bindingCount[interestingCount]
    #stateCount = stateCount[interestingCount]

    
    totalCount = np.array([np.mean(x) for x in componentCount if not math.isnan(np.mean(x))])
    #totalCount = np.array([y for x in componentCount for y in x])
    bindingTotalCount = np.array([y for x in bindingCount for y in x])
    stateTotalCount = np.array([y for x in stateCount for y in x])
    print '----directory: {0}'.format(directory)
    print 'number of models', len(totalCount)
    print 'component summary', np.mean(totalCount), np.std(totalCount), len(totalCount), sum(totalCount)
    print 'binding summary', np.mean(bindingTotalCount), np.std(bindingTotalCount), sum(bindingTotalCount)
    print 'modification summary', np.mean(stateTotalCount), np.std(stateTotalCount), sum(stateTotalCount)
    '''
    plt.clf()
    plt.hist(totalCount, bins=[0, 1, 2, 3, 4, 5, 6, 7])
    plt.xlabel('Number of components', fontsize=20)
    plt.ylabel('Number of molecules', fontsize=20)
    plt.savefig('componentsvsmolecules.png')
    '''
    #plt.clf()
    zeroComponents  = [collections.Counter(x)[0] for x in componentCount]
    
    #print sum([x for x in zeroComponents if x > 7]), len([x for x in zeroComponents if x > 7])
    return totalCount, bindingTotalCount, stateTotalCount, atoarray
    #plt.show()

def getXMLFailures(directory):
    bng = getValidFiles(directory, 'bngl')
    bngxml =  getValidFiles(directory, 'xml')

    bngxml = ['.'.join(x.split('.')[:-1]) + '.bngl' for x in bngxml]

    failures = [x for x in bng if x not in bngxml]
    print len(failures)
    with open(os.path.join(directory, 'failure.dump'), 'wb') as f:
        pickle.dump(failures, f)

def componentDensityPlot():
    '''
    obtains a density plot that compares the distribution of components against three model datasets
    '''

    directory = [('bngTest', 'BNG control set'), ('curated', 'BioModels curated'), ('non_curated', 'BioModels non\n curated')]
    #directory = [('curated', 'BioModels curated')]
    #('new_non_curated', 'BioModels non curated')]
    colors = sns.color_palette("Set1", 3)
    colors = [colors[1], colors[2], colors[0]]
    f, (ax1) = plt.subplots(1, 1, sharex=True, figsize=(6, 3.45))
    f.tight_layout() 
    for color, direct in zip(colors, directory):
        totalCount, bindingCount, modifyCount, atoarray = componentAnalysis(direct[0], 0.1)
        sns.kdeplot(totalCount, color=color, label=direct[1], shade=True, ax=ax1, clip=(0.4, 100), bw=0.2)
        #sns.distplot(bindingCount, color=color, ax=ax2, clip=(-0.1, 8), bw=0.5)
        #sns.distplot(modifyCount, color=color, ax=ax3, clip=(-0.1, 8), bw=0.5)
    plt.xlabel('Number of components', fontsize=22,fontweight='bold')
    #f.text(-0.14,0.5,'Model percentage', fontsize=22,fontweight='bold',va='center', rotation='vertical')
    ax1.set_title('Components/molecule')
    #ax2.set_title('Binding components/molecule')
    #ax3.set_title('Modification components/molecule')
    ax1.set_ylabel('Fraction',fontsize=22,fontweight='bold')
    #ax2.set_ylabel('Fraction',fontsize=22,fontweight='bold')
    #ax3.set_ylabel('Fraction',fontsize=22,fontweight='bold')
    plt.tight_layout()
    ax1.set(xlim=(0,10))
    sns.despine()
    plt.savefig('componentDensity2.pdf',bbox_inches='tight')

    #g = sns.FacetGrid(pandasDistro, row="process", hue="database", margin_titles=True, row_order=processOrder, xlim=(0, 1))
    #g.map(sns.kdeplot, "fraction", shade=True, clip=(0, 1))
    #plt.savefig('componentDensity.pdf',bbox_inches='tight')    



def categorizeStatistics(effort):
    if effort == 0:
        return '0'
    elif effort >=1 and effort <= 5:
        return '1-5'
    elif effort >5 and effort <= 10:
        return '6-10'
    elif effort >10:
        return '>10'


def componentDistroPlot():
    '''
    obtains a bar plot that compares the distribution of components against three model datasets
    '''

    directory = [('bngTest', 'BNG control\n set'), ('curated', 'BioModels\n curated'), ('non_curated', 'BioModels non\n curated')]

    maindb = {'database':[],'category':[],'value':[], 'categorizedScore':[]}
    for direct in directory:
        
        totalCount, bindingCount, modifyCount, atoarray = componentAnalysis(direct[0], 0.1)
        for t,b,m in zip(totalCount,bindingCount, modifyCount):
            maindb['database'].append(direct[1])
            maindb['category'].append('Components')
            maindb['value'].append(t)
            maindb['categorizedScore'].append(categorizeStatistics(t))
            #maindb['database'].append(direct[1])
            #maindb['category'].append('Binding components')
            #maindb['value'].append(b)
            #maindb['database'].append(direct[1])
            #maindb['category'].append('Modification components')
            #maindb['value'].append(m)

    maindb = pandas.DataFrame.from_dict(maindb)
    #g = sns.factorplot(x="database", y="value", row="category", data=maindb, kind="bar",
    #           legend_out=True,aspect=1.9)
    g = sns.factorplot(x="database", y="value", data=maindb, kind="bar",
               legend_out=True,aspect=1.9)

    g.set_xlabels("Database",fontweight='bold',fontsize=28)
    g.set_ylabels("Number",fontweight='bold',fontsize=28)
    g.set_titles("{row_name}",size=24)
    g.despine(left=True)
    g.set_xticklabels(["BNG control\nset","BioModels\ncurated", "BioModels non\n curated"], fontsize=22)
    g.fig.savefig('componentDistro2.pdf',bbox_inches='tight')


def componentvsatomizationPlot():
    directory = [('curated', 'BioModels\n curated'), ('non_curated', 'BioModels non\n curated')]

    maindb = {'database':[],'category':[],'value':[], 'atoscore':[]}
    for direct in directory:
        
        totalCount, bindingCount, modifyCount, atoarray = componentAnalysis(direct[0], 0.1)
        for t,b,m,a in zip(totalCount,bindingCount, modifyCount,atoarray):
            maindb['database'].append(direct[1])
            maindb['category'].append('Components')
            maindb['value'].append(t)
            maindb['atoscore'].append(a)
    maindb = pandas.DataFrame.from_dict(maindb)

    g = sns.FacetGrid(maindb, row="database", hue="database",
        margin_titles=True, xlim=(0, 1), ylim=(0,3.5))
    g.map(sns.kdeplot, "atoscore","value",  shade=True)
    plt.savefig('componentvsato.pdf',bbox_inches='tight')


def processHistogram():
    def getsubplot(axs, index, dimensions):
        if 1 in dimensions:
            return axs[index]
        else:
            return axs[index / dimensions[1]][index%dimensions[0]]
    #directory = [('bnglTest', 'BNG control set'), ('complex2', 'BioModels curated'), ('new_non_curated', 'BioModels non curated')]
    directory = [('bngTest', 'BNG control set'), ('curated', 'BioModels curated'), ('non_curated', 'BioModels non curated')]
    processDistro  = []
    
    cluster = True
    
    pandasDistro = collections.defaultdict(list)
    
    
    
    if not cluster:
        processFile = 'processDistro.dump'
        processOrder = ['Add', 'Delete', 'AddBond', 'DeleteBond', 'StateChange']
    else:
        processFile = 'processDistroCluster.dump'
        processOrder = ['Add/Delete', '(Add/Delete)Bond', 'StateChange']
    '''
    for direct in directory:
        process, pandasc = extractProcessDistribution(direct, cluster=cluster)
        print direct, process
        processDistro.append(process)
        for element in pandasc:
            pandasDistro[element].extend(pandasc[element])
    
    pandasDistro = pandas.DataFrame(data=pandasDistro)
    #print pandasDistro
    with open(processFile, 'wb') as f:
        pickle.dump(processDistro, f)
        pickle.dump(pandasDistro, f)
    '''
    
    with open(processFile, 'rb') as f:
        processDistro = pickle.load(f)
        pandasDistro = pickle.load(f)

    #get only those files that are not entirely syn/del
    tmp = set(pandasDistro.query('process == "Add/Delete" & fraction != 1').file)
    pandasDistro = pandasDistro[pandasDistro.file.isin(tmp)]
    pandasDistro = pandasDistro[pandasDistro.database.isin([x[1] for x in directory])]

    #g = sns.factorplot(x="process", y="fraction", row="database", data=pandasDistro, kind='bar', legend_out=True,order=processOrder,aspect=1.9,palette="Set2")
    g = sns.factorplot(x="database", y="fraction", row="process", data=pandasDistro, kind='bar', legend_out=True,aspect=1.9)

    g.set_xlabels("Process",fontweight='bold',fontsize=28)
    g.set_ylabels("Fraction",fontweight='bold',fontsize=28)
    g.set_titles("{row_name} {row_var}",size=24)
    g.set_xticklabels(fontsize=22)
    g.despine(left=True)

    plt.savefig('processBarChar2.pdf',bbox_inches='tight')

    # place in here the subplot grid size
    
    dimensions = [3, 1]
    f, axs = plt.subplots(dimensions[0], dimensions[1], sharex=True, figsize=(8, 8))
    #sns.set_context('talk', font_scale=1, rc={"lines.linewidth": 2.5})
    colors = sns.color_palette("Set1", 3)
    for color, direct in zip(colors, directory):
        for index, process in enumerate(processOrder):
            legend = False
            if index == len(processOrder) - 1:
                legend = True
            actions = np.array(pandasDistro[pandasDistro.process == process][pandasDistro.database == direct[1]]['fraction'].values)
            sns.kdeplot(actions, shade=True, color=color, label=direct[1], ax=getsubplot(axs, index, dimensions), clip=(0, 1), cumulative=True, legend=legend)
            ax = getsubplot(axs, index, dimensions).set_title(process)
            getsubplot(axs, index, dimensions).set_xlim([0, 1])

    #plt.suptitle('CDF for the fraction of <x> process to the total number of processes in a model for different datasets')
    plt.savefig('processDensity.png',bbox_inches='tight')

    g = sns.FacetGrid(pandasDistro, row="process", hue="database", 
        margin_titles=True, row_order=processOrder, xlim=(0, 1))
    g.map(sns.kdeplot, "fraction", shade=True, clip=(0, 1))
    plt.savefig('processDensityGrid.png',bbox_inches='tight')
    #plt.show()
    #


def rankMoleculeTypes(directory):
    with open(os.path.join(directory, 'moleculeTypeDataSet.dump'), 'rb') as f:
        moleculeTypesArray = pickle.load(f)
    progress = progressbar.ProgressBar()
    moleculeTypeTuples = []
    for element in progress(moleculeTypesArray):
        for molecule in element[0]:
            moleculeTypeTuples.append([molecule.name, len(molecule.components), element[4]])
    moleculeTypesDatabase = pandas.DataFrame(data=moleculeTypeTuples, columns=['molecule', 'size', 'files'])    
    print moleculeTypesDatabase.sort('size', ascending=False).head(30)

def annotationPerAtomizationGroup(directory):
    with open('parsedAnnotations.dump','rb') as f:
        parsedAnnotations = pickle.load(f)

    with open(os.path.join('XMLExamples',directory,'modelAnnotationDictionary.dump'),'rb') as f:
        modelAnnotations = pickle.load(f)

    atoDB = pandas.read_hdf('curatedDB.h5')
    atoDB = atoDB.query('numspecies > 5')
    import pprint
    #low atomization


    lowato = atoDB.query('numspecies < 10 ').index
    lowato = ['XMLExamples/curated/{0}'.format(x.split('/')[-1][0:-4]) for x in lowato]
    #lowatocounterann = collections.Counter([w for x in lowato for y in modelAnnotations[x] for z in modelAnnotations[x][y] for w in parsedAnnotations[z]  if z in parsedAnnotations])
    lowatocounterann = collections.Counter([y for x in lowato for y in modelAnnotations[x] if 'taxonomy' not in y  and 'mamo' not in y])

    print lowatocounterann.most_common(20)

    atoscore = atoDB.query('numspecies >= 10 & numspecies < 30').index
    atoscore = ['XMLExamples/curated/{0}'.format(x.split('/')[-1][0:-4]) for x in atoscore]
    lowatocounterann = collections.Counter([y for x in atoscore for y in modelAnnotations[x]  if 'taxonomy' not in y  and 'mamo' not in y])
    print '---'
    print lowatocounterann.most_common(20)

    atoscore = atoDB.query('numspecies >= 30 & numspecies < 50').index
    atoscore = ['XMLExamples/curated/{0}'.format(x.split('/')[-1][0:-4]) for x in atoscore]

    lowatocounterann = collections.Counter([y for x in atoscore for y in modelAnnotations[x] if 'taxonomy' not in y and 'mamo' not in y])
    print '---'
    print lowatocounterann.most_common(20)


    highato = atoDB.query('numspecies >= 50 & numspecies < 100').index
    highato = ['XMLExamples/curated/{0}'.format(x.split('/')[-1][0:-4]) for x in highato]
    #highatocounterann = collections.Counter([w for x in highato for y in modelAnnotations[x] for z in modelAnnotations[x][y] for w in parsedAnnotations[z]  if z in parsedAnnotations])
    print '+++++', set([x for x in highato if 'http://identifiers.org/go/GO:0019722' in modelAnnotations[x]])

    highatocounterann = collections.Counter([y for x in highato for y in modelAnnotations[x] if 'taxonomy' not in y  and 'mamo' not in y])

    print '---'
    pprint.pprint(highatocounterann.most_common(20))


if __name__ == "__main__":
    #preliminaryAnalysis(directory='curated', directory2='XMLExamples/curated')
    #preliminaryAnalysis(directory='non_curated', directory2='XMLExamples/non_curated')

    #print '---'
    #preliminaryAnalysis(directory='bngTest')
    #rankMoleculeTypes('curated')
    annotationPerAtomizationGroup('curated')
    #processHistogram()

    #componentDensityPlot()    
    #componentvsatomizationPlot()
    #componentDistroPlot()
    '''
    colors = ['r', 'Y', 'b']
    #print processDistro
    fig, ax = plt.subplots()
    #print np.arange(len(processDistro[0].keys()))
    index = 0

    rects = []
    processList = ['StateChange', 'AddBond', 'DeleteBond', 'Add', 'Delete']
    '''
    '''
    for x, y in enumerate(colors):
        rects.append(ax.bar(np.arange(len(processDistro[x].keys()))+0.25*x, [processDistro[x][z][0] for z in processList], 0.25, color=y, yerr=[processDistro[x][z][1] for z in processList]))
        
    

    ax.set_ylabel('Percentage')
    ax.set_title('Process percentage by database')
    ax.set_xticks(np.arange(len(processDistro[x].keys()))+0.35)
    ax.set_xticklabels( [x for x in processList] )

    ax.legend( (rects[0][0], rects[1][0], rects[2][0]), ('BNG control group', 'BioModels curated', 'BioModels non curated') )
    plt.savefig('processBarChar.png')
    '''
    #plt.show()
    #getXMLFailures(directory[0])
    #

    #main('bnglTest', '')
    # rawAnnotationCoverage('XMLExamples/non_curated')
