# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:06:43 2013

@author: proto
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 17:42:31 2011

@author: proto
"""
from pyparsing import Word, Suppress, Optional, alphanums, Group, ZeroOrMore
from copy import deepcopy, copy
import analyzeSBML
import utils.structures as st
from utils.util import logMess
#import biogrid
import marshal
import functools
import utils.pathwaycommons as pwcm
from collections import Counter, defaultdict
import itertools
from atomizerUtils import BindingException

def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = marshal.dumps([args, kwargs])
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


speciesNameGrammar = (Word(alphanums + "_" + ":#-")
                      + Suppress('()') + Optional(Suppress('@' + Word(alphanums + '_-')))) + ZeroOrMore(Suppress('+') + Word(alphanums + "_" + ":#-")
                                                                                                        + Suppress("()") + Optional(Suppress('@' + Word(alphanums + '_-'))))

nameGrammar = Word(alphanums + '_-') + ':'

rateGrammar = Word(alphanums + "()")

grammar = Suppress(Optional(nameGrammar)) + ((Group(speciesNameGrammar) | '0') + Suppress(Optional("<") + "->") +
                                             (Group(speciesNameGrammar) | '0') + Suppress(rateGrammar)) \
    ^ (speciesNameGrammar + Suppress(Optional("<") + "->") + Suppress(rateGrammar))


@memoize
def parseReactions(reaction):

    result = grammar.parseString(reaction).asList()
    if len(result) < 2:
        result = [result, []]
    if '<->' in reaction and len(result[0]) == 1 and len(result[1]) == 2:
        result2 = [result[1], result[0]]
        result = result2
    return result


def resolveDependencyGraph(dependencyGraph, reactant, withModifications=False):
    topCandidate = resolveDependencyGraphHelper(
        dependencyGraph, reactant, [], withModifications)
    return topCandidate


assumptions = defaultdict(set)


def addAssumptions(assumptionType, assumption):
    assumptions[assumptionType].add(assumption)


class CycleError(Exception):

    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, memory):
        self.memory = memory


def getAnnotations(annotation):
    """
    parses a libsbml.XMLAttributes annotation object into a list of annotations
    """
    annotationDictionary = []
    if annotation == [] or annotation is None:
        return []
    for index in range(0, annotation.getNumAttributes()):
        annotationDictionary.append(annotation.getValue(index))
    return annotationDictionary


def getURIFromSBML(moleculeName, parser, filterString=None):
    """
    filters a list of URI's so that we get only filterString  ID's
    """
    annotationList = []
    if parser:
        annotations = parser.getSpeciesAnnotation()
        if annotations[moleculeName]:
            for annotation in annotations[moleculeName]:
                annotationList.extend(getAnnotations(annotation))
    if filterString:
        annotationList = [x for x in annotationList if any(
            filterstr in x for filterstr in filterString)]

    return annotationList


def isInComplexWith(moleculeSet, parser=None):
    """
    given a list of binding candidates, it gets the uniprot ID from annotation and queries
    the pathway commons class to see if there's known binding information for those two
    """
    validPairs = []
    for element in moleculeSet:
        # if element[0] == element[1]:
        #    return []

        name1 = getURIFromSBML(element[0], parser, ['uniprot', 'go'])
        name2 = getURIFromSBML(element[1], parser, ['uniprot', 'go'])
        
        modelAnnotation = parser.extractModelAnnotation()
        modelOrganism = modelAnnotation[
            'BQB_OCCURS_IN'] if 'BQB_OCCURS_IN' in modelAnnotation else None
        molecule1 = name1[0].split('/')[-1] if name1 else element[0]
        molecule2 = name2[0].split('/')[-1] if name2 else element[1]

        simpleOrganism = [x.split('/')[-1]
                          for x in modelOrganism] if modelOrganism else None
        bindingResults = pwcm.queryBioGridByName(
            molecule1, molecule2, simpleOrganism, element[0], element[1])
        if not bindingResults:
            bindingResults = pwcm.queryBioGridByName(
                element[0], element[1], simpleOrganism, None, None)
        #bindingResults = None
        #bindingResults = pwcm.isInComplexWith([element[0], name1], [element[1], name2], organism=modelOrganism)
        if bindingResults:
            validPairs.append(element)
    # use pathway commosn as fallback since its much slower
    if not validPairs:
        for element in moleculeSet:
            # if element[0] == element[1]:
            #    return []

            name1 = getURIFromSBML(element[0], parser, ['uniprot', 'go'])
            name2 = getURIFromSBML(element[1], parser, ['uniprot', 'go'])
            modelAnnotation = parser.extractModelAnnotation()
            modelOrganism = modelAnnotation[
                'BQB_OCCURS_IN'] if 'BQB_OCCURS_IN' in modelAnnotation else None
            bindingResults = pwcm.isInComplexWith(
                [element[0], name1], [element[1], name2], organism=modelOrganism)
            if bindingResults:
                validPairs.append(element)

    return validPairs


@memoize
def resolveDependencyGraphHelper(dependencyGraph, reactant, memory,
                                 withModifications=False):
    result = []
    # if type(reactant) == tuple:
    #    return []
    if reactant not in dependencyGraph or dependencyGraph[reactant] == [] or \
            dependencyGraph[reactant] == [[reactant]]:
        if not withModifications:
            result.append([reactant])
    else:
        for option in dependencyGraph[reactant]:
            tmp = []
            for element in option:
                if element in memory and not withModifications:
                    result.append([element])
                    continue
                elif element in memory:
                    #logMess(
                    #    'ERROR:SCT201', 'dependency cycle detected on {0}'.format(element))
                    raise CycleError(memory)
                baseElement = resolveDependencyGraphHelper(dependencyGraph, element,
                                                           memory + [element], withModifications)
                if baseElement is not None:
                    tmp.extend(baseElement)
            # if not withModifications:
            result.extend(tmp)
            if len(option) == 1 and withModifications and option[0] != reactant:
                result.append((option[0], reactant))
    return result


def addToDependencyGraph(dependencyGraph, label, value):
    if label not in dependencyGraph:
        dependencyGraph[label] = []
    if value not in dependencyGraph[label] and value != []:
        dependencyGraph[label].append(value)


def bindingReactionsAnalysis(dependencyGraph, reaction, classification):
    '''
    adds addBond based reactions based dependencies to the dependency graph
    '''
    totalElements = [item for sublist in reaction for item in sublist]
    for element in totalElements:
        addToDependencyGraph(dependencyGraph, element, [])
        if classification == 'Binding':
            if len(reaction[1]) == 1 and element not in reaction[0]:
                addToDependencyGraph(dependencyGraph, element, reaction[0])
            elif len(reaction[0]) == 1 and element not in reaction[1]:
                addToDependencyGraph(dependencyGraph, element, reaction[1])

def weightDependencyGraph(dependencyGraph):
    @memoize
    def measureGraph(element, path):
        counter = 1
        for x in path:
            if type(x) == list or type(x) == tuple:
                counter += measureGraph(element, x)
            elif x != '0' and x != element:
                counter += 1
        return counter

    weights = []
    for element in dependencyGraph:
        path = resolveDependencyGraph(dependencyGraph, element)
        try:
            path2 = resolveDependencyGraph(dependencyGraph, element, True)
        except CycleError:
            path2 = []
        weight = measureGraph(element, path) + len(path2)
        weights.append([element, weight])

    weights = sorted(weights, key=lambda rule: (rule[1], len(rule[0])))
    return weights


def consolidateDependencyGraph(dependencyGraph, equivalenceTranslator,
                               equivalenceDictionary, sbmlAnalyzer, database, loginformation=True):
    equivalenceTranslator = {}
    def selectBestCandidate(reactant, candidates, dependencyGraph, sbmlAnalyzer,
                            equivalenceTranslator=equivalenceTranslator, equivalenceDictionary=equivalenceDictionary):


        tmpCandidates = []
        modifiedElementsPerCandidate = []
        unevenElements = []
        candidateDict = {}
        for individualAnswer in candidates:
            try:
                tmpAnswer = []
                flag = True
                if len(individualAnswer) == 1 and individualAnswer[0] == reactant:
                    continue
                modifiedElements = []
                for chemical in individualAnswer:

                    # we cannot handle tuple naming conventions for now
                    if type(chemical) == tuple:
                        flag = False
                        continue
                    # associate elements in the candidate description with their
                    # modified version
                    rootChemical = resolveDependencyGraph(
                        dependencyGraph, chemical)
                    mod = resolveDependencyGraph(dependencyGraph, chemical, True)
                    if mod != []:
                        modifiedElements.extend(mod)
                    for element in rootChemical:
                        if len(element) == 1 and type(element[0]) == tuple:
                            continue
                        if element == chemical:
                            tmpAnswer.append(chemical)
                        elif type(element) == tuple:
                            tmpAnswer.append(element)
                        else:
                            tmpAnswer.append(element[0])
                modifiedElementsPerCandidate.append(modifiedElements)
                if flag:
                    tmpAnswer = sorted(tmpAnswer)
                    tmpCandidates.append(tmpAnswer)
            except CycleError:
                if loginformation:
                    logMess('ERROR:SCT221', '{0}:Dependency cycle found when mapping molecule to candidate {1}'.format(reactant, individualAnswer))
                continue
        # we cannot handle tuple naming conventions for now
        if len(tmpCandidates) == 0:
            # logMess('CRITICAL:Atomization', 'I dont know how to process these candidates and I have no \
            # way to make an educated guess. Politely refusing to translate
            # {0}={1}.'.format(reactant, candidates))
            return None, None, None
        originalTmpCandidates = deepcopy(tmpCandidates)
        # if we have more than one modified element for a single reactant
        # we can try to  choose the one that is most similar to the original
        # reactant
        # FIXME:Fails if there is a double modification
        newModifiedElements = {}
        #modifiedElementsCounter = Counter()
        modifiedElementsCounters = [Counter() for x in range(len(candidates))]
        # keep track of how many times we need to modify elements in the candidate description
        # FIXME: This only keeps track of the stuff in the fist candidates list
        for idx, modifiedElementsInCandidate in enumerate(modifiedElementsPerCandidate):
            for element in modifiedElementsInCandidate:
                if element[0] not in newModifiedElements or element[1] == reactant:
                    newModifiedElements[element[0]] = element[1]
                modifiedElementsCounters[idx][element[0]] += 1

        # actually modify elements and store final version in tmpCandidates
        # if tmpCandidates[1:] == tmpCandidates[:-1] or len(tmpCandidates) ==
        # 1:

        for tmpCandidate, modifiedElementsCounter in zip(tmpCandidates, modifiedElementsCounters):
            flag = True
            while flag:
                flag = False
                for idx, chemical in enumerate(tmpCandidate):
                    if modifiedElementsCounter[chemical] > 0:
                        modifiedElementsCounter[chemical] -= 1
                        tmpCandidate[idx] = newModifiedElements[chemical]
                        flag = True
                        break
        candidateDict = {
            tuple(x): y for x, y in zip(tmpCandidates, candidates)}
        bcan = []
        btmp = []
        borig = []
        # filter out those dependencies to the 0 element

        # if this is related to the zero element
        if len(tmpCandidates) == 1 and tmpCandidates[0] == ['0']:
            return ['0'], None, None

        for candidate, tmpcandidate, originaltmpcandidate in zip(candidates, tmpCandidates, originalTmpCandidates):
            if originaltmpcandidate != ['0']:
                bcan.append(candidate)
                btmp.append(tmpcandidate)
                borig.append(originaltmpcandidate)
        candidates = bcan
        tmpCandidates = btmp
        originalTmpCandidates = borig

        if len(tmpCandidates) == 0:
            return None, None, None

        # FIXME: I have no idea wtf this is doing so im commenting it out. i
        # think it's old code that is no longer ncessary
        '''
        # update candidate chemical references to their modified version if required
        if len(tmpCandidates) > 1:
            # temporal solution for defaulting to the first alternative
            totalElements = [y for x in tmpCandidates for y in x]
            elementDict = {}
            for word in totalElements:
                if word not in elementDict:
                    elementDict[word] = 0
                elementDict[word] += 1
            newTmpCandidates = [[]]
            for element in elementDict:
                if elementDict[element] % len(tmpCandidates) == 0:
                    newTmpCandidates[0].append(element)
                #elif elementDict[element] % len(tmpCandidates) != 0 and re.search('(_|^){0}(_|$)'.format(element),reactant):
                #    newTmpCandidates[0].append(element)
                #    unevenElements.append([element])
                else:
                    logMess('WARNING:Atomization', 'Are these actually the same? {0}={1}.'.format(reactant,candidates))
                    unevenElements.append(element)
            flag = True
            # FIXME:this should be done on newtmpCandidates instead of tmpcandidates
            while flag:
                flag = False
                for idx, chemical in enumerate(tmpCandidates[0]):
                    if chemical in newModifiedElements: #and newModifiedElements[chemical] in reactant:
                        tmpCandidates[0][idx] = newModifiedElements[chemical]
                        flag = True
                        break
        '''
        # if all the candidates are about modification changes to a complex
        # then try to do it through lexical analysis
        if all([len(candidate) == 1 for candidate in candidates]) and \
                candidates[0][0] != reactant and len(tmpCandidates[0]) > 1:
            if reactant is not None:
                pass

            # analyze based on standard modifications
            #lexCandidate, translationKeys, tmpequivalenceTranslator = sbmlAnalyzer.analyzeSpeciesModification(candidates[0][0], reactant, originalTmpCandidates[0])
            # print '++++'
            lexCandidate, translationKeys, tmpequivalenceTranslator = sbmlAnalyzer.analyzeSpeciesModification2(
                candidates[0][0], reactant, originalTmpCandidates[0])
            # lexCandidate, translationKeys, tmpequivalenceTranslator = sbmlAnalyzer.analyzeSpeciesModification(candidates[0][0], reactant, tmpCandidates[0])            # FIXME: this is iffy. is it always an append modification? could be prepend
            #lexCandidate = None
            if lexCandidate is not None:
                lexCandidate = tmpCandidates[0][
                    originalTmpCandidates[0].index(lexCandidate)]
                lexCandidateModification = lexCandidate + translationKeys[0]
                for element in tmpequivalenceTranslator:
                    if element not in equivalenceTranslator:
                        equivalenceTranslator[element] = []
                equivalenceTranslator[element].append(
                    (lexCandidate, lexCandidateModification))
                dependencyGraph[lexCandidateModification] = [[lexCandidate]]
                while lexCandidate in tmpCandidates[0]:
                    tmpCandidates[0].remove(lexCandidate)
                    tmpCandidates[0].append(lexCandidateModification)
                    break
                return [tmpCandidates[0]], unevenElements, candidates

            else:
                fuzzyCandidateMatch = None
                '''
                if nothing else works and we know the result is a bimolecular
                complex and we know which are the basic reactants then try to
                do fuzzy string matching between the two.
                TODO: extend this to more than 2 molecule complexes.
                '''
                if len(tmpCandidates[0]) == 2:
                    tmpmolecules = []
                    tmpmolecules.extend(originalTmpCandidates[0])
                    tmpmolecules.extend(tmpCandidates[0])
                    # FIXME: Fuzzy artificial reaction is using old methods. Try to fix this
                    # or maybe not, no one was using it and when it was used it was wrong
                    # fuzzyCandidateMatch = sbmlAnalyzer.fuzzyArtificialReaction(originalTmpCandidates[0],[reactant],tmpmolecules)
                    fuzzyCandidateMatch = None
                if fuzzyCandidateMatch is not None:
                    # logMess('INFO:Atomization', 'Used fuzzy string matching from {0} to {1}'.format(reactant, fuzzyCandidateMatch))
                    return [fuzzyCandidateMatch], unevenElements, candidates
                else:
                    # map based on greedy matching
                    greedyMatch = sbmlAnalyzer.greedyModificationMatching(
                        reactant, dependencyGraph.keys())
                    if greedyMatch not in [-1, -2]:
                        return selectBestCandidate(reactant, [greedyMatch], dependencyGraph, sbmlAnalyzer)[0], unevenElements, candidates

                    # last ditch attempt using straighforward lexical analysis
                    tmpDependency, tmpEquivalence = sbmlAnalyzer.findClosestModification(
                        [reactant], dependencyGraph.keys(), database)
                    if reactant in tmpDependency and tmpDependency[reactant] in tmpCandidates[0]:
                        for element in tmpDependency:
                            if element not in dependencyGraph:
                                dependencyGraph[
                                    element] = tmpDependency[element]
                        for element in tmpEquivalence:
                            if element not in equivalenceDictionary:
                                equivalenceDictionary[element] = []
                            for equivalence in tmpEquivalence[element]:
                                if equivalence[0] not in equivalenceDictionary[element]:
                                    equivalenceDictionary[
                                        element].append(equivalence[0])
                        if len(tmpDependency.keys()) > 0:
                            return tmpDependency[reactant], unevenElements, candidates
                    # XXX: be careful of this change. This basically forces changes to happen
                    # the ive no idea whats going on branch
                    #modificationCandidates = {}
                    # if modificationCandidates == {}:
                    if len(tmpCandidates) != 1:
                        if not database.softConstraints:
                            if loginformation:
                                logMess('ERROR:SCT213', '{0}:Atomizer needs user information to determine which element is being modified among components {1}={2}.'.format(
                                reactant, candidates, tmpCandidates))
                            # print database.userLabelDictionary
                            return None, None, None
                    else:
                        if not database.softConstraints:
                            if loginformation:
                                logMess('ERROR:SCT212', '{0}:Atomizer needs user information to determine which element is being modified among component species {1}={2}.'.format(
                                reactant, candidates, tmpCandidates))
                            return None, None, None
                        # print database.userLabelDictionary')
                    # return [tmpCandidates[0]], unevenElements

        elif len(tmpCandidates) > 1:
            # all candidates are equal/consistent
            if all(sorted(x) == sorted(tmpCandidates[0]) for x in tmpCandidates):
                tmpCandidates = [tmpCandidates[0]]
            elif reactant in database.alternativeDependencyGraph and loginformation:
                # candidates contradict each other but we have naming convention information in alternativeDependencyGraph
                if not all(sorted(x) == sorted(originalTmpCandidates[0]) for x in originalTmpCandidates):
                    if loginformation:
                        logMess('INFO:SCT001', '{0}:Using lexical analysis since stoichiometry gives non-consistent information naming({1})!=stoichiometry({2})'.format(reactant,
                                                                                                                                                           database.alternativeDependencyGraph[reactant][0], tmpCandidates))

                # else:
                #    print database.alternativeDependencyGraph[reactant],tmpCandidates,reactant
                #    logMess('INFO:Atomization', 'Using lexical analysis for species {0} =  {1} since stoichiometry gave conflicting information {2}'.format(reactant,
                # database.alternativeDependencyGraph[reactant][0],
                # tmpCandidates))

                # fallback to naming conventions
                candidate = database.alternativeDependencyGraph[reactant]
                # resolve naming convention candidate to its basic components
                # (molecule types)
                namingTmpCandidates = selectBestCandidate(
                    reactant, [candidate[0]], dependencyGraph, sbmlAnalyzer)[0]
                if not namingTmpCandidates:
                    logMess('ERROR:SCT211', '{0}:Cannot converge to solution, conflicting definitions {1}={2}'.format(
                            reactant, tmpCandidates, originalTmpCandidates))                    
                    return None, None, None
                if not any([sorted(subcandidate) == sorted(namingTmpCandidates[0]) for subcandidate in tmpCandidates]):
                    if loginformation:
                        logMess('WARNING:SCT112', '{0}:Stoichiometry analysis result in non self-consistent definitions but conflicts with lexical analysis stoichiometry({1})!= naming({2}). Selecting lexical analysis'.format(reactant,
                                                                                                                                                                                                                                          tmpCandidates, namingTmpCandidates))
                    addAssumptions('lexicalVsstoch', (reactant, ('lexical', str(
                        namingTmpCandidates)), ('stoch', str(tmpCandidates)), ('original', str(originalTmpCandidates))))

                tmpCandidates = namingTmpCandidates
                if loginformation:
                    database.alternativeDependencyGraph[reactant] = tmpCandidates
            elif all(sorted(x) == sorted(originalTmpCandidates[0]) for x in originalTmpCandidates):
                #the basic elements are the same but we are having trouble matching modifciations together
                sortedCandidates = sorted([([y for y in x if y in reactant], i) for i, x in enumerate(
                    tmpCandidates)], key=lambda z: [len(z[0]), sum([len(w) for w in z[0]])], reverse=True)
                if loginformation:
                    logMess('WARNING:SCT113', '{0}:candidates {1} agree on the basic components but naming conventions cannot determine  specific modifications. Selecting {2} based on longest partial match'.format(
                        reactant, tmpCandidates, tmpCandidates[sortedCandidates[0][1]]))
                replacementCandidate = [tmpCandidates[sortedCandidates[0][1]]]
                addAssumptions('lexicalVsstoch', (reactant, ('current', str(replacementCandidate)), ('alternatives', str(
                    [x for x in tmpCandidates if x != replacementCandidate[0]])), ('original', str(originalTmpCandidates))))
                tmpCandidates = replacementCandidate
            else:
                tmpCandidates2 = [x for x in tmpCandidates if all(y not in x for y in database.constructedSpecies)]
                # if we had constructed species disregard those since they are introducing noise
                if len(tmpCandidates2) > 0 and len(tmpCandidates) != len(tmpCandidates2):
                    return selectBestCandidate(reactant, tmpCandidates2, dependencyGraph, sbmlAnalyzer)
                elif len(tmpCandidates2) == 0:
                    #the differences is between species that we created so its the LAE fault. Just choose one.
                    tmpCandidates.sort(key=len)
                    tmpCandidates = [tmpCandidates[0]]
                else:
                    if loginformation:
                        logMess('ERROR:SCT211', '{0}:Cannot converge to solution, conflicting definitions {1}={2}'.format(
                        reactant, tmpCandidates, originalTmpCandidates))
                    return None, None, None
        elif reactant in database.alternativeDependencyGraph and loginformation:
            # there is one stoichionetry candidate but the naming convention
            # and the stoichionetry dotn agree
            if tmpCandidates[0] != database.alternativeDependencyGraph[reactant][0]:
                # make sure the naming convention is resolved to basic
                # omponents
                candidate = database.alternativeDependencyGraph[reactant]
                # this is to avoid recursion
                if loginformation:
                    del database.alternativeDependencyGraph[reactant]
                namingtmpCandidates = selectBestCandidate(
                    reactant, [candidate[0]], dependencyGraph, sbmlAnalyzer)[0]

                # if they still disagree print error and use stoichiometry
                if namingtmpCandidates and tmpCandidates[0] != namingtmpCandidates[0]:

                    if loginformation:
                        if namingtmpCandidates[0][0] in database.constructedSpecies:
                            namingTmpCandidates = tmpCandidates

                        else:
                            database.alternativeDependencyGraph[reactant] = namingtmpCandidates
                            logMess('WARNING:SCT111', '{0}:conflicting definitions between stoichiometry ({1}) and naming conventions {2}. Choosing the latter'.format(
                                    reactant, tmpCandidates[0], database.alternativeDependencyGraph[reactant]))
                    tmpCandidates = namingtmpCandidates
                    addAssumptions('lexicalVsstoch', (reactant, ('stoch', str(tmpCandidates)), ('lexical', str(
                        namingtmpCandidates)), ('original', str(originalTmpCandidates))))
                    for element in tmpCandidates[0]:
                        if element not in prunnedDependencyGraph:
                            # elemental species that were not used anywhere
                            # else but for those entries discovered through
                            # naming conventions
                            prunnedDependencyGraph[element] = []
                elif not namingtmpCandidates:
                    if loginformation:
                        logMess('WARNING:SCT121','{0}:could not resolve naming({1}) into a viable compositional candidate. choosing stoichiometry({2})'.format(reactant,candidate,tmpCandidates[0]))
        originalCandidateName = candidateDict[tuple(tmpCandidates[0])] if tuple(
            tmpCandidates[0]) in candidateDict else None
        return [tmpCandidates[0]], unevenElements, originalCandidateName

    prunnedDependencyGraph = deepcopy(dependencyGraph)

    tempMergedDependencyGraph = deepcopy(prunnedDependencyGraph)
    for element in database.alternativeDependencyGraph:
        if element in tempMergedDependencyGraph:
            tempMergedDependencyGraph[element].extend(database.alternativeDependencyGraph[element])
    weights = weightDependencyGraph(tempMergedDependencyGraph)



    #raise Exception

    unevenElementDict = {}
    for element in weights:
        candidates = [x for x in prunnedDependencyGraph[element[0]]]
        if len(candidates) == 1 and type(candidates[0][0]) == tuple:
            prunnedDependencyGraph[element[0]] = []
        if len(candidates) >= 1:
            candidates, uneven, originalCandidate = selectBestCandidate(
                element[0], candidates, prunnedDependencyGraph, sbmlAnalyzer)

            #except CycleError:
            #    candidates = None
            #    uneven = []
            if uneven != []:
                unevenElementDict[element[0]] = (uneven)
        if candidates is None:
            prunnedDependencyGraph[element[0]] = []
        else:
            prunnedDependencyGraph[element[0]] = [
                sorted(x) for x in candidates]

    weights = weightDependencyGraph(prunnedDependencyGraph)
    return prunnedDependencyGraph, weights, unevenElementDict, equivalenceTranslator


def identifyReaction(equivalenceDictionary, originalElement, modifiedElement):
    for classification in equivalenceDictionary:
        if set([originalElement, modifiedElement]) in [set(x) for x in equivalenceDictionary[classification]]:
            return classification
    return None


def createEmptySpecies(name):
    species = st.Species()
    molecule = st.Molecule(name)
    species.addMolecule(molecule)
    return species


def addStateToComponent(species, moleculeName, componentName, state):
    for molecule in species.molecules:
        if moleculeName == molecule.name:
            for component in molecule.components:
                if componentName == component.name:
                    tmp = component.activeState
                    if state not in component.states:
                        component.addState(state)
                    elif state in component.states:
                        component.setActiveState(state)
                    return tmp


def addComponentToMolecule(species, moleculeName, componentName):
    for molecule in species.molecules:
        if moleculeName == molecule.name:
            if componentName not in [x.name for x in molecule.components]:
                component = st.Component(componentName)
                molecule.addComponent(component)
                return True
    return False


def addBondToComponent(species, moleculeName, componentName, bond, priority=1):
    order = 1
    for molecule in species.molecules:
        if moleculeName == molecule.name:
            if priority == order or len([x for x in species.molecules if x.name == moleculeName]) == order:
                for component in molecule.components:
                    if componentName == component.name:
                        # if we are adding a second bond to the same component
                        # it actually means that we have two components with the
                        # same name
                        if len(component.bonds) == 0:
                            component.addBond(bond)
                            return
                        else:
                            newComponent = st.Component(componentName)
                            newComponent.addBond(bond)
                            molecule.addComponent(newComponent)
                            return
            else:
                order += 1


def solveComplexBinding(totalComplex, pathwaycommonsFlag, parser, compositionEntry):
    '''
    given two binding complexes it will attempt to find the ways in which they bind using different criteria

    '''
    def sortMolecules(array, reverse):
        return sorted(array, key=lambda molecule: (len(molecule.components),
                                                   len([
                                                       x for x in molecule.components if x.activeState not in [0, '0']]),
                                                   len(str(molecule)), str(molecule)), reverse=reverse)

    def getBiggestMolecule(array):
        sortedMolecule = sortMolecules(array, reverse=False)

        #sortedMolecule = sorted(sortedMolecule, key=lambda rule: len(rule.components))

        return sortedMolecule[-1]

    def getNamedMolecule(array, name):

        for molecule in sortMolecules(array, True):
            if molecule.name == name:
                return molecule
            elif molecule.trueName == name:
                return molecule

    names1 = [str(x.trueName) for x in totalComplex[0] ]
    names2 = [str(x.trueName) for x in totalComplex[1] ]
    bioGridDict = {}
    # find all pairs of molecules
    comb = set([tuple(sorted([x, y])) for x in names1 for y in names2])
    comb2 = set([tuple(sorted([x,y])) for x in compositionEntry for y in compositionEntry])

    dbPair = set([])
    combTemp = set()

    # search pathway commons for binding candidates
    if pathwaycommonsFlag:
        dbPair = isInComplexWith(comb, parser)
    else:
        for element in comb:
            if element[0].upper() in bioGridDict and element[1] in bioGridDict[element[0].upper()] or \
                    element[1].upper() in bioGridDict and element[0] in bioGridDict[element[1].upper()]:
                #logMess('INFO:ATO001', 'Biogrid info: {0}:{1}'.format(element[0], element[1]))
                dbPair.add((element[0], element[1]))
        # elif pathwaycommonsFlag:
        #    if pwcm.isInComplexWith(element[0], element[1]):
        #        dbPair.add((element[0], element[1]))
    dbPair = list(dbPair)

    if dbPair != []:
        mol1 = mol2 = None
        # select the best candidate if there's many ways to bind (in general
        # one that doesn't overlap with an already exising pair)
        finalDBpair = []
        if len(dbPair) > 1:

            for element in dbPair:
                mset1 = Counter(element)
                mset2 = Counter(names1)
                mset3 = Counter(names2)
                intersection1 = mset1 & mset2
                intersection2 = mset1 & mset3
                intersection1 = list(intersection1.elements())
                intersection2 = list(intersection2.elements())
                if len(intersection1) < 2 and len(intersection2) < 2:
                    finalDBpair.append(element)
        if len(finalDBpair) > 0:
            dbPair = finalDBpair

        if len(dbPair) > 1:

            # @FIXME: getNamedMolecule should never receive parameters that cause it to return null, but somehow that's what is happening
            # when you receive a malformed user definition file. The error
            # should be caught way before we reach this point
            tmpComplexSubset1 = [getNamedMolecule(totalComplex[0], element[
                                                  0]) for element in dbPair if getNamedMolecule(totalComplex[0], element[0]) is not None]
            if not tmpComplexSubset1:
                tmpComplexSubset1 = [getNamedMolecule(totalComplex[0], element[
                                                      1]) for element in dbPair if getNamedMolecule(totalComplex[0], element[1]) is not None]
                tmpComplexSubset2 = [getNamedMolecule(totalComplex[1], element[
                                                      0]) for element in dbPair if getNamedMolecule(totalComplex[1], element[0]) is not None]
            else:
                tmpComplexSubset2 = [getNamedMolecule(totalComplex[1], element[
                                                      1]) for element in dbPair if getNamedMolecule(totalComplex[1], element[1]) is not None]

            mol1 = getBiggestMolecule(tmpComplexSubset1)
            mol2 = getBiggestMolecule(tmpComplexSubset2)
            logMess('INFO:ATO002', "According to BioGrid/Pathwaycommons there's more than one way to bind {0} and {1} together: {2}. Defaulting to {3}-{4}".format(names1, names2, dbPair,
                                                                                                                                                                   mol1.name, mol2.name))

        else:
            mol1 = getNamedMolecule(totalComplex[0], dbPair[0][0])
            if not mol1:
                mol1 = getNamedMolecule(totalComplex[1], dbPair[0][0])
                mol2 = getNamedMolecule(totalComplex[0], dbPair[0][1])

            else:
                mol2 = getNamedMolecule(totalComplex[1], dbPair[0][1])
                if not mol2:
                    mol1 = getNamedMolecule(totalComplex[1], dbPair[0][0])
                    mol2 = getNamedMolecule(totalComplex[0], dbPair[0][1])

            logMess('INFO:ATO001', 'Binding information found in BioGrid/Pathwaycommons for for {0}-{1}'.format(mol1.name, mol2.name))

    else:

        #mol1 = getBiggestMolecule(totalComplex[0])
        #mol2 = getBiggestMolecule(totalComplex[1])
        '''
        if pathwaycommonsFlag:
            logMess('ERROR:ATO201', "We don't know how {0} and {1} bind together and there's no relevant BioGrid/Pathway commons information. Not atomizing".format(
                [x.name for x in totalComplex[0]], [x.name for x in totalComplex[1]]))
            # addAssumptions('unknownBond',(mol1.name,mol2.name))
        else:

            logMess('ERROR:ATO202', "We don't know how {0} and {1} bind together. Not atomizing".format(
                [x.name for x in totalComplex[0]], [x.name for x in totalComplex[1]]))
            # addAssumptions('unknownBond',(mol1.name,mol2.name))
        '''
        raise BindingException(
            '{0}-{1}'.format(sorted([x.name for x in totalComplex[0]]), sorted([x.name for x in totalComplex[1]])), comb)

    return mol1, mol2


def getComplexationComponents2(moleculeName, species, bioGridFlag, pathwaycommonsFlag=False,
                               parser=None, bondSeeding=[], bondExclusion=[], database=None):
    '''
    method used during the atomization process. It determines how molecules
    in a species bind together
    '''

    def sortMolecules(array, reverse):
        return sorted(array, key=lambda molecule: (len(molecule.components),
                                                   len([
                                                       x for x in molecule.components if x.activeState not in [0, '0']]),
                                                   len(str(molecule)), str(molecule)),
                      reverse=reverse)

    def getBiggestMolecule(array):
        sortedMolecule = sortMolecules(array, reverse=False)

        #sortedMolecule = sorted(sortedMolecule, key=lambda rule: len(rule.components))

        return sortedMolecule[-1]

    def getNamedMolecule(array, name):

        for molecule in sortMolecules(array, True):
            if molecule.name == name:
                return molecule
    speciesDict = {}
    # this array will contain all molecules that bind together
    pairedMolecules = []
    for x in sortMolecules(species.molecules, reverse=True):
        for y in x.components:
            if y.name not in speciesDict:
                speciesDict[y.name] = []
            speciesDict[y.name].append(x)
    # this array wil contain all molecules that dont bind to anything
    orphanedMolecules = [x for x in species.molecules]
    # seed the list of pairs from the seeds
    pairedMolecules = copy(bondSeeding)
    if bondSeeding:
        orphanedMolecules = [
            x for x in orphanedMolecules for y in bondSeeding if x not in y]

    # determine how molecules bind together
    redundantBonds = []

    for x in sortMolecules(species.molecules, reverse=True):
        for component in [y for y in x.components if y.name.lower()
                          in speciesDict.keys()]:
            if x.name.lower() in speciesDict:
                if(x in speciesDict[component.name.lower()]) and component.name in [y.name.lower() for y
                                                                                    in speciesDict[x.name.lower()]]:
                    for mol in speciesDict[x.name.lower()]:
                        if mol.name.lower() == component.name and x != mol and x in \
                                speciesDict[component.name]:
                            speciesDict[x.name.lower()].remove(mol)
                            speciesDict[component.name].remove(x)
                            if x not in orphanedMolecules and mol not in orphanedMolecules:
                                # FIXME: is it necessary to remove double bonds
                                # in complexes?
                                
                                lhs = set([])
                                rhs = set([])
                                repeatedFlag = False
                                for pair in pairedMolecules:
                                    

                                    if x in pair:
                                        lhs.add(pair[0])
                                        lhs.add(pair[1])
                                    elif mol in pair:
                                        rhs.add(pair[0])
                                        rhs.add(pair[1])
                                    # is this particular pair of molecules bound together?
                                    if x in pair and mol in pair:
                                        repeatedFlag = True
                                        break
                                # this pair already exists
                                if repeatedFlag:
                                    continue
                                redundantBonds.append([x, mol])
                                intersection = lhs.intersection(rhs)
                                redundantBonds[-1].extend(list(intersection))
                                if len(redundantBonds[-1]) < 3:
                                    redundantBonds.pop()
                                # continue
                            if [x, mol] not in bondSeeding and [mol, x] not in bondSeeding and [x, mol] not in bondExclusion and [mol, x] not in bondExclusion:
                                pairedMolecules.append([x, mol])
                            if x in orphanedMolecules:
                                orphanedMolecules.remove(x)
                            if mol in orphanedMolecules:
                                orphanedMolecules.remove(mol)

    if len(redundantBonds) > 0:
        for x in redundantBonds:
            addAssumptions(
                'redundantBonds', tuple(sorted([y.name for y in x])))
            addAssumptions(
                'redundantBondsMolecules', (tuple(sorted([y.name for y in x])), moleculeName))
            logMess('WARNING:CTX001', 'Redundant bonds detected between molecules {0} in species {1}'.format(
                [y.name for y in x], moleculeName))
    totalComplex = [set(x) for x in pairedMolecules]
    isContinuousFlag = True

    # iterate over orphaned and find unidirectional interactions
    # e.g. if a molecule has a previous known interaction with the
    # same kind of molecule, even if it has no available components
    # e.g. k-mers`

    for element in speciesDict:
        for individualMolecule in speciesDict[element]:
            if individualMolecule in orphanedMolecules:
                candidatePartner = [x for x in species.molecules if x.name.lower(
                ) == element and x != individualMolecule]
                if len(candidatePartner) == 1:
                    pairedMolecules.append(
                        [candidatePartner[0], individualMolecule])
                    orphanedMolecules.remove(individualMolecule)
    # determine which pairs form a continuous chain

    while isContinuousFlag:
        isContinuousFlag = False
        for idx in range(0, len(totalComplex) - 1):
            for idx2 in range(idx + 1, len(totalComplex)):
                if len([x for x in totalComplex[idx] if x in totalComplex[idx2]]) > 0:
                    totalComplex[idx] = totalComplex[
                        idx].union(totalComplex[idx2])
                    totalComplex.pop(idx2)
                    isContinuousFlag = True
                    break
            if isContinuousFlag:
                break
    # now we process those molecules where we need to create a new component
    for element in orphanedMolecules:

        for mol1 in species.molecules:
            # when adding orphaned molecules make sure it's not already in
            # the list
            if mol1 == element and mol1 not in set().union(*totalComplex):
                totalComplex.append(set([mol1]))
    # now we process for those molecules we are not sure how do they bind
    while len(totalComplex) > 1:

        if len(totalComplex[0]) == 1 and len(totalComplex[1]) == 1:
            mol1 = list(totalComplex[0])[0]
            mol2 = list(totalComplex[1])[0]
        else:
            mol1, mol2 = solveComplexBinding(totalComplex, pathwaycommonsFlag, parser, database.prunnedDependencyGraph[moleculeName][0])
        pairedMolecules.append([mol1, mol2])
        totalComplex[0] = totalComplex[0].union(totalComplex[1])
        totalComplex.pop(1)
    # totalComplex.extend(orphanedMolecules)
    return pairedMolecules


def getTrueTag(dependencyGraph, molecule):
    '''
    given any modified or basic element it returns its basic
    name
    '''
    if dependencyGraph[molecule] == []:
        return molecule
    elif dependencyGraph[molecule][0][0] == molecule:
        return molecule
    else:
        return getTrueTag(dependencyGraph, dependencyGraph[molecule][0][0])


def createCatalysisRBM(dependencyGraph, element, translator, reactionProperties,
                       equivalenceDictionary, sbmlAnalyzer, database):
    '''
    if it's a catalysis reaction create a new component/state
    '''
    if dependencyGraph[element[0]][0][0] == element[0]:
        if element[0] not in translator:
            translator[element[0]] = createEmptySpecies(element[0])
    else:
        componentStateArray = []
        tmp = element[0]
        existingComponents = []
        memory = []
        forceActivationSwitch = False
        while dependencyGraph[tmp] != []:
            # what kind of catalysis are we dealing with
            # classification = identifyReaction(
            #                                  equivalenceDictionary,
            #                                  dependencyGraph[tmp][0][0],tmp)

            classifications = None
            if not classifications:
                classifications = identifyReaction(
                    equivalenceDictionary,
                    dependencyGraph[tmp][0][0], tmp)
                classifications = classifications if classifications in reactionProperties else None
                if classifications is not None:
                    classifications = [classifications]
            if not classifications:
                classifications = sbmlAnalyzer.findMatchingModification(
                    tmp, dependencyGraph[tmp][0][0])

            if not classifications:
                classifications = sbmlAnalyzer.findMatchingModification(
                    element[0], dependencyGraph[tmp][0][0])

            # if we know what classification it is then add the corresponding
            # components and states
            if classifications is not None:
                for classification in classifications:
                    componentStateArray.append(
                        reactionProperties[classification])
                    # classificationArray.append([classification,
                    #                            tmp,dependencyGraph[tmp]
                    #                            [0][0]])
                    existingComponents.append(
                        reactionProperties[classification][0])
            # if we don't know we can create a force 1:1 modification
            elif database.forceModificationFlag and classifications is None and not forceActivationSwitch:
                forceActivationSwitch = True
                baseName = getTrueTag(dependencyGraph,
                                      dependencyGraph[element[0]][0][0])

                species = createEmptySpecies(baseName)

                componentStateArray.append(['{0}'.format(tmp), tmp])
                logMess('WARNING:LAE002', 'adding forced transformation: {0}:{1}:{2}'.format(
                    baseName, dependencyGraph[element[0]], element[0]))
                # return
            # bail out if we couldn't figure out what modification it is
            elif classifications is None:
                logMess('DEBUG:MSC001', 'unregistered modification: {0}:{1}'.format(
                    element[0], dependencyGraph[element[0]]))
            memory.append(tmp)
            tmp = dependencyGraph[tmp][0][0]
            if tmp in memory:
                raise CycleError(memory)
        baseName = getTrueTag(
            dependencyGraph, dependencyGraph[element[0]][0][0])

        species = createEmptySpecies(baseName)
        # use the already existing structure if its in the
        # translator,otherwise empty
        if baseName in translator:
            species = translator[baseName]
        #modifiedSpecies = deepcopy(translator[dependencyGraph[element[0]][0][0]])

        # modified species needs to start from the base speceis sine componentStateArray should contain the full set of modifications
        # check that this works correctly for double modifications
        modifiedSpecies = deepcopy(translator[baseName])
        # this counter is here for multi level modification events (e.g. double
        # phosporylation)
        modificationCounter = {
            componentState[0]: 2 for componentState in componentStateArray}
        for componentState in componentStateArray:
            addComponentToMolecule(species, baseName, componentState[0])
            addComponentToMolecule(
                modifiedSpecies, baseName, componentState[0])
            tmp = addStateToComponent(species, baseName,
                                      componentState[0], componentState[1])
            if tmp == componentState[1]:
                addStateToComponent(species, baseName, componentState[0],
                                    componentState[1] + componentState[1])
            # this modification was already activated so create a second
            # modification component
            if addStateToComponent(modifiedSpecies, baseName,
                                   componentState[0], componentState[1]) == componentState[1]:
                componentName = '{0}{1}'.format(
                    componentState[0],  modificationCounter[componentState[0]])
                modificationCounter[componentState[0]] += 1
                addComponentToMolecule(
                    modifiedSpecies, baseName, componentName)
                addStateToComponent(modifiedSpecies, baseName, componentName,
                                    componentState[1])
            addStateToComponent(species, baseName, componentState[0], '0')
        # update the base species
        if len(componentStateArray) > 0:
            translator[baseName] = deepcopy(species)
            translator[element[0]] = modifiedSpecies


globalNumberGenerator = []


def getBondNumber(molecule1, molecule2):
    """
    keeps a model-level registry of of all the molecule pairs and returns a unique index
    """
    moleculeList = tuple(sorted([molecule1, molecule2]))
    if moleculeList not in globalNumberGenerator:
        globalNumberGenerator.append(moleculeList)

    return globalNumberGenerator.index(moleculeList)


def createBindingRBM(element, translator, dependencyGraph, bioGridFlag, pathwaycommonsFlag, parser, database):
    species = st.Species()

    # go over the sct and reuse existing stuff

    for molecule in dependencyGraph[element[0]][0]:
        if molecule in translator:
            tmpSpecies = translator[molecule]
            if molecule != getTrueTag(dependencyGraph, molecule):
                original = translator[getTrueTag(dependencyGraph, molecule)]
                updateSpecies(tmpSpecies, original.molecules[0])
            if tmpSpecies.molecules[0].name in database.constructedSpecies:
                tmpSpecies.molecules[0].trueName = molecule
            else:
                tmpSpecies.molecules[0].trueName = tmpSpecies.molecules[0].name
            species.addMolecule(deepcopy(tmpSpecies.molecules[0]))
        else:
            mol = st.Molecule(molecule)
            mol.trueName = molecule
            #dependencyGraph[molecule] = deepcopy(mol)
            species.addMolecule(mol)
    dependencyGraphCounter = Counter(dependencyGraph[element[0]][0])

    # XXX: this wont work for species with more than one molecule with the
    # same name
    changeFlag = False
    partialBonds = defaultdict(list)
    for partialUserEntry in database.partialUserLabelDictionary:
        partialCounter = Counter(partialUserEntry)
        if all([partialCounter[x] <= dependencyGraphCounter[x] for x in partialCounter]):
            changeFlag = True
            for molecule in database.partialUserLabelDictionary[partialUserEntry].molecules:
                for molecule2 in species.molecules:
                    if molecule.name == molecule2.name:
                        for component in molecule.components:
                            for bond in component.bonds:
                                if molecule2.name not in [x.name for x in partialBonds[bond]]:
                                    partialBonds[bond].append(molecule2)
                        '''
                        for component in molecule.components:
                            component2 = [x for x in molecule2.components if x.name == component.name]
                            # component already exists in species template
                            if component2:
                                if component.bonds:
                                    component2[0].bonds = component.bonds
                            else:
                                molecule2.addComponent(deepcopy(component))
                        '''

    bondSeeding = [partialBonds[x] for x in partialBonds if x > 0]
    bondExclusion = [partialBonds[x] for x in partialBonds if x < 0]
    # how do things bind together?
    moleculePairsList = getComplexationComponents2(
        element[0], species, bioGridFlag, pathwaycommonsFlag, parser, bondSeeding, bondExclusion, database)
    #moleculeCount = Counter([y for x in moleculePairsList for y in x])
    # print moleculeCount
    #moleculePairsList = [sorted(x) for x in moleculePairsList]
    #moleculePairsList.sort(key=lambda x: [-moleculeCount[x[0]],(str(x[0]), x[0],str(x[1]),x[1])])
    # TODO: update basic molecules with new components
    # translator[molecule[0].name].molecules[0].components.append(deepcopy(newComponent1))
    # translator[molecule[1].name].molecules[0].components.append(deepcopy(newComponent2))
    moleculeCounter = defaultdict(list)
    for molecule in moleculePairsList:
        flag = False

        # create an index on m0 and m1 depending on their name and repeats in
        # the species
        if molecule[0] not in moleculeCounter[molecule[0].name]:
            moleculeCounter[molecule[0].name].append(molecule[0])
        if molecule[1] not in moleculeCounter[molecule[1].name]:
            moleculeCounter[molecule[1].name].append(molecule[1])
        m0index = moleculeCounter[molecule[0].name].index(molecule[0])
        m1index = moleculeCounter[molecule[1].name].index(molecule[1])
        bondIdx = getBondNumber('{0}{1}'.format(
            molecule[0].name, m0index), '{0}{1}'.format(molecule[1].name, m1index))
        # add bonds where binding components already exist and they are not
        # occupied
        for component in molecule[0].components:
            if component.name == molecule[1].name.lower() and \
                    len(component.bonds) == 0:
                component.bonds.append(bondIdx)
                flag = True
                break
        if not flag:
            # create components if they dont exist already.
            # Add a bond afterwards
            newComponent1 = st.Component(molecule[1].name.lower())

            molecule[0].components.append(newComponent1)

            if newComponent1.name not in [x.name for x in translator[molecule[0].name].molecules[0].
                                          components]:
                translator[molecule[0].name].molecules[
                    0].components.append(deepcopy(newComponent1))

            molecule[0].components[-1].bonds.append(bondIdx)
        flag = False
        # same thing for the other member of the bond
        for component in molecule[1].components:
            if component.name == molecule[0].name.lower() and len(component.bonds) == 0:
                component.bonds.append(bondIdx)
                flag = True
                break
        if not flag:
            newComponent2 = st.Component(molecule[0].name.lower())
            molecule[1].components.append(newComponent2)
            if molecule[0].name != molecule[1].name:
                if newComponent2.name not in [x.name for x in translator[molecule[1].name].molecules[0].
                                              components]:
                    translator[molecule[1].name].molecules[
                        0].components.append(deepcopy(newComponent2))
            molecule[1].components[-1].bonds.append(bondIdx)

    # update the translator
    translator[element[0]] = species


def atomize(dependencyGraph, weights, translator, reactionProperties,
            equivalenceDictionary, bioGridFlag, sbmlAnalyzer, database, parser):
    '''
    The atomizer's main methods. Receives a dependency graph
    '''
    redrawflag = True
    loops = 0
    while redrawflag and loops < 10:
        loops +=1
        bindingCounter = Counter()
        bindingFailureDict = {}

        for idx, element in enumerate(weights):
            # 0 molecule
            if element[0] == '0':
                continue
            # user defined molecules to be the zero molecule
            if dependencyGraph[element[0]] == [['0']]:
                zeroSpecies = st.Species()
                zeroMolecule = st.Molecule('0')
                zeroSpecies.addMolecule(zeroMolecule)
                translator[element[0]] = zeroSpecies
                continue
            # undivisible molecules

            elif dependencyGraph[element[0]] == []:
                if element[0] not in translator:
                    translator[element[0]] = createEmptySpecies(element[0])
            else:
                if len(dependencyGraph[element[0]][0]) == 1:
                    # catalysis
                    createCatalysisRBM(dependencyGraph, element, translator, reactionProperties,
                                       equivalenceDictionary, sbmlAnalyzer, database)
                else:
                    try:
                        createBindingRBM(element, translator, dependencyGraph,
                                         bioGridFlag, database.pathwaycommons, parser, database)

                    except BindingException as e:
                        for c in e.combinations:
                            bindingCounter[c] += 1
                            bindingFailureDict[element[0]] = e.combinations
                        logMess('DEBUG:ATO003', "We don't know how {0} binds together in complex {1}. Not atomizing".format(
                                e.value, element[0]))
            
                        # there awas an issue during binding, don't atomize
                        translator[element[0]] = createEmptySpecies(element[0])
        
        # evaluate species that weren't bound properly and see if we can get information from all over the model to find the right binding partner
        bindingTroubleLog = defaultdict(list)
        modifiedPairs = set()

        redrawflag = False
        for molecule in bindingFailureDict:
            bindingWinner = defaultdict(list)
            for candidateTuple in bindingFailureDict[molecule]:
                bindingWinner[bindingCounter[candidateTuple]].append(candidateTuple)
            bestBindingCandidates = bindingWinner[max(bindingWinner.keys())]
            if len(bestBindingCandidates) > 1:
                bindingTroubleLog[tuple(sorted(bestBindingCandidates))].append(molecule)
            else:
                bindingPair = bestBindingCandidates[0]
                if bindingPair not in modifiedPairs:
                    modifiedPairs.add(bindingPair)
                else:
                    continue
                c1 = st.Component(bindingPair[1].lower())
                c2 = st.Component(bindingPair[0].lower())
                molecule1 = translator[translator[bindingPair[0]].molecules[0].name].molecules[0]
                molecule2 = translator[translator[bindingPair[1]].molecules[0].name].molecules[0]

                molecule1.addComponent(c1)
                molecule2.addComponent(c2)
                redrawflag = True
                logMess('INFO:ATO031','Determining that {0} binds together based on frequency of the bond in the reaction network.'.format(bindingPair))
    for trouble in bindingTroubleLog:
        logMess('ERROR:ATO202','{0}:We need information to resolve the bond structure of these complexes . \
Please choose among the possible binding candidates that had the most observed frequency in the reaction network or provide a new one:{1}'.format(bindingTroubleLog[trouble],trouble))

    

def updateSpecies(species, referenceMolecule):
    flag = False
    for moleculeStructure in species.molecules:
        if moleculeStructure.name == referenceMolecule.name:
            for component in referenceMolecule.components:
                count = [x.name for x in referenceMolecule.components].count(
                    component.name)
                count -= [x.name for x in moleculeStructure.components].count(
                    component.name)
                newComponent = st.Component(component.name)
                if len(component.states) > 0:
                    newComponent.addState('0')
                if count > 0:
                    for _ in range(0, count):
                        #just make a copy of the reference component and set active state to 0
                        moleculeStructure.addComponent(deepcopy(component))
                        moleculeStructure.components[-1].setActiveState('0')
                elif count < 0:
                    for _ in range(0, -count):
                        #FIXME: does not fully copy the states
                        referenceMolecule.addComponent(deepcopy(newComponent))
                        flag = True
                elif count == 0:
                    localComponents = [
                        x for x in moleculeStructure.components if x.name == component.name]
                    referenceComponents = [
                        x for x in referenceMolecule.components if x.name == component.name]
                    if [x.states for x in localComponents] != [x.states for x in referenceComponents]:
                        for lc in localComponents:
                            for rc in referenceComponents:
                                for ls in lc.states:
                                    if ls not in rc.states:
                                        rc.addState(ls, update=False)
                                for rs in rc.states:
                                    if rs not in lc.states:
                                        lc.addState(rs, update=False)

            for component in moleculeStructure.components:
                count = [x.name for x in referenceMolecule.components].count(
                    component.name)
                count -= [x.name for x in moleculeStructure.components].count(
                    component.name)
                newComponent = st.Component(component.name)
                if len(component.states) > 0:
                    newComponent.addState(component.states[0])
                    newComponent.addState('0')
                if count > 0:
                    for idx in range(0, count):
                        moleculeStructure.addComponent(deepcopy(newComponent))
                elif count < 0:
                    for idx in range(0, -count):
                        referenceMolecule.addComponent(deepcopy(newComponent))
                        flag = True

    return flag


def propagateChanges(translator, dependencyGraph):
    flag = True
    while flag:
        flag = False
        for dependency in dependencyGraph:
            if dependencyGraph[dependency] == []:
                continue
            for molecule in dependencyGraph[dependency][0]:
                try:
                    if updateSpecies(translator[dependency],
                                     translator[getTrueTag(dependencyGraph,
                                                           molecule)].molecules[0]):
                        flag = True
                except:
                    logMess('CRITICAL:Program',
                            'Species is not being properly propagated')
                    flag = False


def fillSCTwithAnnotationInformation(orphanedSpecies, annotationDict, database, logResults=True,tentativeFlag=False):
        # annotation handling
    exactMatches = defaultdict(list)
    partialMatches = defaultdict(list)
    strongIntersectionMatches = defaultdict(list)
    intersectionMatches = defaultdict(list)
    # iterate over all pairs of orphaned species
    for combinationParticle in itertools.combinations(orphanedSpecies, 2):
            # compare annotations
        if combinationParticle[0] in annotationDict and combinationParticle[1] in annotationDict:

            sortedPair = sorted(list(combinationParticle), key=len)
            # get unary keys
            unaryAnnotation1 = [y for x in annotationDict[combinationParticle[0]] for y in annotationDict[
                combinationParticle[0]][x] if x in ['BQM_IS_DESCRIBED_BY', 'BQB_IS_VERSION_OF', 'BQB_IS','BQB_ENCODES'] and ('uniprot' in y or 'chebi' in y)]
            unaryAnnotation2 = [y for x in annotationDict[combinationParticle[1]] for y in annotationDict[
                combinationParticle[1]][x] if x in ['BQM_IS_DESCRIBED_BY', 'BQB_IS_VERSION_OF', 'BQB_IS','BQB_ENCODES'] and ('uniprot' in y or 'chebi' in y)]

            # get compositional keys
            compositionalAnnotation1 = [y for x in annotationDict[combinationParticle[0]] for y in annotationDict[
                combinationParticle[0]][x] if x in ['BQB_HAS_PART', 'BQB_HAS_VERSION'] and ('uniprot' in y or 'chebi' in y)]
            compositionalAnnotation2 = [y for x in annotationDict[combinationParticle[1]] for y in annotationDict[
                combinationParticle[1]][x] if x in ['BQB_HAS_PART', 'BQB_HAS_VERSION'] and ('uniprot' in y or 'chebi' in y)]
            # unary keys match
            if any([x in unaryAnnotation2 for x in unaryAnnotation1]):

                exactMatches[sortedPair[1]].append([sortedPair[0]])
            # one composes the other
            elif any([x in compositionalAnnotation1 for x in unaryAnnotation2]):
                #if combinationParticle[0] not in partialMatches:
                #    partialMatches[combinationParticle[0]].append([])
                partialMatches[combinationParticle[0]].append([combinationParticle[1]])

            elif any([x in compositionalAnnotation2 for x in unaryAnnotation1]):
                #if combinationParticle[1] not in partialMatches:
                #    partialMatches[combinationParticle[1]].append([])
                partialMatches[combinationParticle[1]].append([combinationParticle[0]])
            elif set(compositionalAnnotation1) == set(compositionalAnnotation2) and len([x in compositionalAnnotation2 for x in compositionalAnnotation1]) > 0:
                strongIntersectionMatches[sortedPair[1]].append([sortedPair[0]])
            # they intersect
            elif any([x in compositionalAnnotation2 for x in compositionalAnnotation1]):

                intersectionMatches[sortedPair[1]].append([sortedPair[0]])
                # intersectionMatches[combinationParticle[0]].append(combinationParticle[1])
    # create unary groups


    exactMatches = consolidateDependencyGraph(
        dict(exactMatches), {}, {}, database.sbmlAnalyzer, database, loginformation=False)[0]

    if logResults:
        for x in [y for y in exactMatches if len(exactMatches[y]) > 0]:
            if not tentativeFlag:
                logMess('INFO:ANN001', '{0}: can be the same as {1} according to annotation information. No action was taken'.format(
                    x, exactMatches[x]))
            else:

                if not (x in database.dependencyGraph and exactMatches[x][0][0] in database.dependencyGraph and database.dependencyGraph[x] == database.dependencyGraph[exactMatches[x][0][0]]):
                    logMess('WARNING:ANN101', '{0}: was determined to be the same as {1} according to annotation information. Please confirm from user information'.format(
                    x, exactMatches[x]))

    # create strong intersection groups

    strongIntersectionMatches = {x: strongIntersectionMatches[x] for x in strongIntersectionMatches if x not in partialMatches}
    strongIntersectionMatches.update(exactMatches)
    strongIntersectionMatches = consolidateDependencyGraph(dict(strongIntersectionMatches), {}, {}, database.sbmlAnalyzer, database, loginformation=False)[0]
    if logResults:
        for x in [y for y in strongIntersectionMatches if len(strongIntersectionMatches[y]) > 0]:
            if x not in exactMatches:
                if not tentativeFlag:
                    logMess('INFO:ANN002', '{0}: can exactly match {1} according to annotation information. No action was taken'.format(
                        x, strongIntersectionMatches[x]))
                else:
                    if not(x in database.dependencyGraph and strongIntersectionMatches[x][0][0] in database.dependencyGraph and \
                        database.dependencyGraph[x] == database.dependencyGraph[strongIntersectionMatches[x][0][0]]):
                        logMess('WARNING:ANN101', '{0}: was determined to exactly match {1} according to annotation information. Please confirm from user information'.format(
                        x, strongIntersectionMatches[x]))
    # create partial intersection groups
    intersectionMatches = {x: intersectionMatches[x] for x in intersectionMatches if x not in partialMatches and x not in strongIntersectionMatches}
    intersectionMatches.update(exactMatches)


    intersectionMatches = consolidateDependencyGraph(dict(intersectionMatches), {}, {}, database.sbmlAnalyzer, database, loginformation=False)[0]
    if logResults:
        for x in intersectionMatches:
            if x not in exactMatches:
                logMess('INFO:ANN002', '{0}: was determined to be partially match {1} according to annotation information.'.format(
                    x, intersectionMatches[x]))

    partialMatches = consolidateDependencyGraph(
        dict(partialMatches), {}, {}, database.sbmlAnalyzer, database, loginformation=False)[0]

    if logResults:
        for x in partialMatches:
            if partialMatches[x] != []:
                logMess('INFO:ANN003', '{0}: is thought to be partially composed of reactants {1} according to annotation information. Please verify stoichiometry'.format(
                    x, partialMatches[x]))

    # validAnnotationPairs.sort()

    intersectionMatches.update(strongIntersectionMatches)

    return intersectionMatches, partialMatches


def createSpeciesCompositionGraph(parser, database, configurationFile, namingConventions,
                                  speciesEquivalences=None, bioGridFlag=False):
    '''
    Creates the species composition table (SCT).

    It first does stoichiometry analysis, then lexical...
    '''
    _, rules, _ = parser.getReactions(atomize=True)
    molecules, _, _, _, _, _ = parser.getSpecies()
    database.sbmlAnalyzer = \
        analyzeSBML.SBMLAnalyzer(
            parser, configurationFile, namingConventions, speciesEquivalences, conservationOfMass=True)

    # classify reactions
    database.classifications, equivalenceTranslator, database.eequivalenceTranslator,\
        indirectEquivalenceTranslator, \
        adhocLabelDictionary, lexicalDependencyGraph, userEquivalenceTranslator = database.sbmlAnalyzer.classifyReactions(
            rules, molecules, {})
    database.reactionProperties = database.sbmlAnalyzer.getReactionProperties()
    # user defined and lexical analysis naming conventions are stored here
    database.reactionProperties.update(adhocLabelDictionary)

    database.translator, database.userLabelDictionary, \
        database.lexicalLabelDictionary, database.partialUserLabelDictionary = database.sbmlAnalyzer.getUserDefinedComplexes()
    database.dependencyGraph = {}
    database.alternativeDependencyGraph = {}

    # fill in the annotation dictionary
    database.annotationDict = parser.getFullAnnotation()
    #database.annotationDict = {}

    # ###dependency graph
    # binding reactions
    for reaction, classification in zip(rules, database.classifications):
        bindingReactionsAnalysis(database.dependencyGraph,
                                 list(parseReactions(reaction)), classification)
    # lexical dependency graph contains lexically induced binding compositions. atomizer gives preference to binding obtained this way as opposed to stoichiometry
    # stronger bounds on stoichiometry based binding can be defined in
    # reactionDefinitions.json.
    for element in lexicalDependencyGraph:

        if element in database.dependencyGraph and element not in database.userLabelDictionary:
            if len(lexicalDependencyGraph[element]) == 0:
                continue
            '''
            oldDependency = database.dependencyGraph[element]
            
            if sorted(lexicalDependencyGraph[element][0]) in [sorted(x) for x in oldDependency]:
                # if len(oldDependency) > 1:
                #    logMess('DEBUG:Atomization', 'Species {0} was confirmed to be {1} based on lexical information'.format(element,lexicalDependencyGraph[element]))
                database.dependencyGraph[
                    element] = lexicalDependencyGraph[element]
            else:
                # logMess('INFO:Atomization', 'Species {0} was determined to be {1} instead of {2} based on \
                # lexical information'.format(element,
                # lexicalDependencyGraph[element], oldDependency))
            '''
            if database.dependencyGraph[element] != []:
                database.alternativeDependencyGraph[element] = lexicalDependencyGraph[element]
            else:
                logMess('INFO:LAE009','{0}: being set to be a modification of cosntructed species {1}'.format(element,lexicalDependencyGraph[element][0]))
                addToDependencyGraph(database.dependencyGraph,element,lexicalDependencyGraph[element][0])
        else:
            database.constructedSpecies.add(element)
            database.dependencyGraph[element] = lexicalDependencyGraph[element]
        # Check if I'm using a molecule that hasn't been used yet
        for dependencyCandidate in database.dependencyGraph[element]:
            for molecule in [x for x in dependencyCandidate if x not in database.dependencyGraph]:
                # this is a species that was not originally in the model. in case theres conflict later this is 
                # to indicate it is given less priority
                database.dependencyGraph[molecule] = []

    # user defined transformations
    for key in userEquivalenceTranslator:
        for namingEquivalence in userEquivalenceTranslator[key]:
            baseElement = min(namingEquivalence, key=len)
            modElement = max(namingEquivalence, key=len)
            if baseElement not in database.dependencyGraph:
                database.dependencyGraph[baseElement] = []
            addToDependencyGraph(database.dependencyGraph, modElement,
                                 [baseElement])


    # database.eequivalence translator contains 1:1 equivalences
    # FIXME: do we need this update step or is it enough with the later one?
    # catalysis reactions
    '''
    for key in database.eequivalenceTranslator:
        for namingEquivalence in database.eequivalenceTranslator[key]:
            baseElement = min(namingEquivalence, key=len)
            modElement = max(namingEquivalence, key=len)
            if key != 'Binding':
                if baseElement not in database.dependencyGraph or database.dependencyGraph[baseElement] == []:
                    if modElement not in database.dependencyGraph or database.dependencyGraph[modElement] == []:
                        database.dependencyGraph[baseElement] = []
                    # do we have a meaningful reverse dependence?
                    # elif all([baseElement not in x for x in database.dependencyGraph[modElement]]):
                    #    addToDependencyGraph(database.dependencyGraph,baseElement,[modElement])
                    #    continue

                        if baseElement in database.annotationDict and modElement in database.annotationDict:
                            baseSet = set([y for x in database.annotationDict[
                                          baseElement] for y in database.annotationDict[baseElement][x]])
                            modSet = set([y for x in database.annotationDict[
                                         modElement] for y in database.annotationDict[modElement][x]])
                            if len(baseSet.intersection(modSet)) > 0 or len(baseSet) == 0 or len(modSet) == 0:
                                addToDependencyGraph(database.dependencyGraph, modElement,
                                                     [baseElement])
                            else:
                                logMess("ERROR:ANN201", "{0} and {1} have a direct correspondence according to reaction information however their annotations are completely different.".format(
                                    baseElement, modElement))
                        else:
                            addToDependencyGraph(database.dependencyGraph, modElement,
                                                 [baseElement])
    '''
    # include user label information.
    for element in database.userLabelDictionary:
        if database.userLabelDictionary[element] in [0, [(0,)]]:
            database.dependencyGraph[element] = ['0']
        elif len(database.userLabelDictionary[element][0]) == 0 or element == \
                database.userLabelDictionary[element][0][0]:

            database.dependencyGraph[element] = []
        else:
            database.dependencyGraph[element] = [list(
                database.userLabelDictionary[element][0])]
            # If the user is introducing a new molecule term, add it to the SCT
            if database.userLabelDictionary[element][0][0] not in database.dependencyGraph:
                database.dependencyGraph[
                    database.userLabelDictionary[element][0][0]] = []

    # add species elements defined by the user into the naming convention
    # definition
    molecules.extend(['{0}()'.format(
        x) for x in database.userLabelDictionary if '{0}()'.format(x) not in molecules])
    # recalculate 1:1 equivalences now with binding information
    _, _, database.eequivalenceTranslator2,\
        _, adhocLabelDictionary, _, _ = database.sbmlAnalyzer.classifyReactions(
            rules, molecules, database.dependencyGraph)
    database.reactionProperties.update(adhocLabelDictionary)

    # update catalysis equivalences
    # catalysis reactions
    for key in database.eequivalenceTranslator2:
        for namingEquivalence in database.eequivalenceTranslator2[key]:

            baseElement = min(namingEquivalence, key=len)
            modElement = max(namingEquivalence, key=len)
            # dont overwrite user information
            if key != 'Binding' and modElement not in database.userLabelDictionary:
                if baseElement not in database.dependencyGraph:
                    database.constructedSpecies.add(baseElement)
                    database.dependencyGraph[baseElement] = []
                if modElement not in database.dependencyGraph or not [True for x in database.dependencyGraph[modElement] if baseElement in x and len(x) > 1]:
                    if baseElement in database.annotationDict and modElement in database.annotationDict:
                        baseSet = set([y for x in database.annotationDict[baseElement] for y in database.annotationDict[baseElement][x]])
                        modSet = set([y for x in database.annotationDict[
                                      modElement] for y in database.annotationDict[modElement][x]])
                        if len(baseSet.intersection(modSet)) > 0 or len(baseSet) == 0 or len(modSet) == 0:
                            if modElement not in database.dependencyGraph:
                                # if the entry doesnt exist from previous information accept this
                                addToDependencyGraph(database.dependencyGraph, modElement, [baseElement])
                            else:
                                # otherwise add it to the lexical repository
                                addToDependencyGraph(database.alternativeDependencyGraph, modElement, [baseElement])
                        else:
                            baseDB = set([x.split('/')[-2] for x in baseSet if 'identifiers.org' in x])
                            modDB = set([x.split('/')[-2] for x in modSet if 'identifiers.org' in x])
                            # it is still ok if they each refer to different databases
                            if len(baseDB.intersection(modDB)) == 0:
                                if modElement not in database.dependencyGraph:
                                    # if the entry doesnt exist from previous information accept this
                                    addToDependencyGraph(database.dependencyGraph, modElement, [baseElement])
                                else:
                                    # otherwise add it to the lexical repository
                                    addToDependencyGraph(database.alternativeDependencyGraph, modElement, [baseElement])
                            else:
                                logMess("WARNING:ANN201", "{0} and {1} have a direct correspondence according to reaction \
information however their annotations are completely different.".format(baseElement, modElement))
                    else:
                        addToDependencyGraph(database.dependencyGraph, modElement,
                                             [baseElement])
                else:
                    logMess('WARNING:ATO114', 'Definition conflict between binding information {0} and lexical analyis {1} for molecule {2},\
choosing binding'.format(database.dependencyGraph[modElement], baseElement, modElement))
    # non lexical-analysis catalysis reactions
    if database.forceModificationFlag:
        for reaction, classification in zip(rules, database.classifications):
            preaction = list(parseReactions(reaction))
            if len(preaction[0]) == 1 and len(preaction[1]) == 1:
                if (preaction[0][0] in [0, '0']) or (preaction[1][0] in [0, '0']):
                    continue
                if preaction[1][0].lower() in preaction[0][0].lower() or len(preaction[1][0]) < len(preaction[0][0]):
                    base = preaction[1][0]
                    mod = preaction[0][0]
                else:
                    mod = preaction[1][0]
                    base = preaction[0][0]
                if database.dependencyGraph[mod] == [] and mod not in database.userLabelDictionary:

                    if base in database.userLabelDictionary and \
                            database.userLabelDictionary[base] == 0:
                        continue
                    if mod in database.userLabelDictionary and \
                            database.userLabelDictionary[mod] == 0:
                        continue
                    if [mod] in database.dependencyGraph[base]:
                        continue

                    # can we just match it up through existing species instead of forcing a modification?
                    greedyMatch = database.sbmlAnalyzer.greedyModificationMatching(
                        mod, database.dependencyGraph.keys())

                    if greedyMatch not in [-1,-2, []]:
                        database.dependencyGraph[mod] = [greedyMatch]
                        if mod in database.alternativeDependencyGraph:
                            del database.alternativeDependencyGraph[mod]
                        logMess('INFO:LAE006','{0}: Mapped to {1} using lexical analysis/greedy matching'.format(mod, greedyMatch))
                        continue
                    

                    # if the annotations have no overlap whatsoever don't force
                    # this modifications
                    if base in database.annotationDict and mod in database.annotationDict:
                        baseSet = set([y for x in database.annotationDict[
                                      base] for y in database.annotationDict[base][x]])
                        modSet = set(
                            [y for x in database.annotationDict[mod] for y in database.annotationDict[mod][x]])
                        if(len(baseSet.intersection(modSet))) == 0 and len(baseSet) > 0 and len(modSet) > 0:
                            baseDB = set([x.split('/')[-2] for x in baseSet if 'identifiers.org' in x])
                            modDB = set([x.split('/')[-2] for x in modSet if 'identifiers.org' in x])
                            #we stil ahve to check that they both reference the same database
                            if len(baseDB.intersection(modDB)) > 0:
                                logMess("WARNING:ANN201", "{0} and {1} have a direct correspondence according to reaction \
information however their annotations are completely different.".format(base, mod))
                                continue
                    database.dependencyGraph[mod] = [[base]]

    '''
    #complex catalysis reactions
    for key in indirectEquivalenceTranslator:
        #first remove these entries from the dependencyGraph since
        #they are not true bindingReactions
        for namingEquivalence in indirectEquivalenceTranslator[key]:
            removedElement = ''
            tmp3 = deepcopy(namingEquivalence[1])
            if tmp3 in database.dependencyGraph[namingEquivalence[0][0]]:
                removedElement = namingEquivalence[0][0]
            elif tmp3 in database.dependencyGraph[namingEquivalence[0][1]]:
                removedElement = namingEquivalence[0][1]

            else:
                tmp3.reverse()
                if tmp3 in database.dependencyGraph[namingEquivalence[0][0]]:
                    removedElement = namingEquivalence[0][0]

                elif tmp3 in database.dependencyGraph[namingEquivalence[0][1]]:
                    removedElement = namingEquivalence[0][1]


            #then add the new, true dependencies
            #if its not supposed to be a basic element
            tmp = [x for x in namingEquivalence[1] if x not in namingEquivalence[2]]
            tmp.extend([x for x in namingEquivalence[2] if x not in namingEquivalence[1]])
            tmp2 = deepcopy(tmp)
            tmp2.reverse()
            ##TODO: map back for the elements in namingEquivalence[2]
            if tmp not in database.dependencyGraph[namingEquivalence[3][0]] \
                and tmp2 not in database.dependencyGraph[namingEquivalence[3][0]]:
                if sorted(tmp) == sorted(tmp3):
                    continue
                if all(x in database.dependencyGraph for x in tmp):
                    if removedElement in database.dependencyGraph:
                        database.dependencyGraph[removedElement].remove(tmp3)
                    logMess('INFO:Atomization','Removing {0}={1} and adding {2}={3} instead\
 from the dependency list since we determined it is not a true binding reaction based on lexical analysis'\
                    .format(removedElement,tmp3,namingEquivalence[3][0],tmp))
                    database.dependencyGraph[namingEquivalence[3][0]] = [tmp]
                else:
                    logMess('WARNING:Atomization','We determined that {0}={1} based on lexical analysis instead of \
{2}={3} (stoichiometry) but one of the constituent components in {1} is not a molecule so no action was taken'.format(namingEquivalence[3][0],
tmp,removedElement,tmp3))
    #user defined stuff
'''

    # stuff obtained from string similarity analysis
    for element in database.lexicalLabelDictionary:
        # similarity analysis has less priority than anything we discovered
        # before
        if element in database.dependencyGraph and \
                len(database.dependencyGraph[element]) > 0:
            continue

        if len(database.lexicalLabelDictionary[element][0]) == 0 or element == \
                database.lexicalLabelDictionary[element][0][0]:
            database.constructedSpecies.add(element)
            addToDependencyGraph(database.dependencyGraph, element, [])
        else:
            # logMess('INFO:Atomization', ' added induced speciesStructure {0}={1}'
            #         .format(element, database.lexicalLabelDictionary[element][0]))
            database.dependencyGraph[element] = [list(
                database.lexicalLabelDictionary[element][0])]



    # Now let's go for annotation analysis and last resort stuff on the remaining orphaned molecules

    strippedMolecules = [x.strip('()') for x in molecules]
    orphanedSpecies = [
        x for x in strippedMolecules if x not in database.dependencyGraph or database.dependencyGraph[x] == []]
    orphanedSpecies.extend([x for x in database.dependencyGraph if database.dependencyGraph[
                           x] == [] and x not in orphanedSpecies])

    # Fill SCT with annotations for those species that still dont have any
    # mapping

    annotationDependencyGraph, _ = fillSCTwithAnnotationInformation(
        orphanedSpecies, database.annotationDict, database)

    # use an empty dictionary if we wish to turn off annotation information in atomization
    #annotationDependencyGraph = {}

    for annotatedSpecies in annotationDependencyGraph:
        if len(annotationDependencyGraph[annotatedSpecies]) > 0 and annotatedSpecies not in database.userLabelDictionary:
            addToDependencyGraph(
                database.dependencyGraph, annotatedSpecies, annotationDependencyGraph[annotatedSpecies][0])
            logMess('INFO:ANN004', 'Added equivalence from annotation information {0}={1}'.format(annotatedSpecies,
                                                                                                  annotationDependencyGraph[annotatedSpecies][0]))
            for element in annotationDependencyGraph[annotatedSpecies][0]:
                # in case one of the compositional elements is not yet in the
                # dependency graph
                if element not in database.dependencyGraph:
                    addToDependencyGraph(database.dependencyGraph, element, [])

    nonOrphanedSpecies = [x for x in strippedMolecules if x not in orphanedSpecies]


    annotationDependencyGraph, _ = fillSCTwithAnnotationInformation(
        nonOrphanedSpecies, database.annotationDict, database,tentativeFlag=True)


    orphanedSpecies = [
        x for x in strippedMolecules if x not in database.dependencyGraph or database.dependencyGraph[x] == []]
    orphanedSpecies.extend([x for x in database.dependencyGraph if database.dependencyGraph[
                           x] == [] and x not in orphanedSpecies])


    orphanedSpecies.extend(database.constructedSpecies)
    strippedMolecules.extend([x.strip('()') for x in database.constructedSpecies])
    # TODO: merge both lists and use them as a tiebreaker for consolidation
    #completeAnnotationDependencyGraph, completePartialMatches = fillSCTwithAnnotationInformation(strippedMolecules, annotationDict, database, False)
    # pure lexical analysis for the remaining orphaned molecules

    tmpDependency, database.tmpEquivalence = database.sbmlAnalyzer.findClosestModification(
        orphanedSpecies, strippedMolecules, database)


    for species in tmpDependency:
        if species not in database.userLabelDictionary:
            if tmpDependency[species] == []:
                addToDependencyGraph(database.dependencyGraph, species, [])
            for instance in tmpDependency[species]:
                addToDependencyGraph(database.dependencyGraph, species, instance)
                if len(instance) == 1 and instance[0] not in database.dependencyGraph:
                    addToDependencyGraph(database.dependencyGraph, instance[0], [])



    orphanedSpecies = [
        x for x in strippedMolecules if x not in database.dependencyGraph or database.dependencyGraph[x] == []]
    orphanedSpecies.extend([x for x in database.dependencyGraph if database.dependencyGraph[
                           x] == [] and x not in orphanedSpecies])
    orphanedSpecies.extend(database.constructedSpecies)

    # greedy lexical analysis for the remaining orhpaned species
    for reactant in orphanedSpecies:
        greedyMatch = database.sbmlAnalyzer.greedyModificationMatching(
            reactant, database.dependencyGraph.keys())
        if greedyMatch not in [-1,-2, []]:
            database.dependencyGraph[reactant] = [greedyMatch]
            logMess('INFO:LAE006','Mapped {0} to {1} using lexical analysis/greedy matching'.format(reactant, greedyMatch))
    if len(database.constructedSpecies) > 0:
        logMess('WARNING:SCT131','The following species names do not appear in the original model but where created to have more appropiate naming conventions: [{0}]'.format(','.join(database.constructedSpecies)))
    # initialize and remove zero elements

    database.prunnedDependencyGraph, database.weights, unevenElementDict, database.artificialEquivalenceTranslator = \
        consolidateDependencyGraph(database.dependencyGraph, equivalenceTranslator,
                                   database.eequivalenceTranslator, database.sbmlAnalyzer, database)

    return database


def sanityCheck(translator):
    '''
    checks for critical atomization errors like isomorphism
    '''
    stringrep = {x: str(translator[x]) for x in translator}
    repeats = set()
    for key in range(0, len(translator.keys()) - 1):
        for key2 in range(key + 1, len(translator.keys())):
            if stringrep[translator.keys()[key]] == stringrep[translator.keys()[key2]]:
                repeats.add((translator.keys()[key], translator.keys()[key2]))
    for repeat in repeats:
        logMess('ERROR:MSC201', 'Elements {0} and {1} produce\
            the same translation. Emptying {1}.'.format(repeat[0], repeat[1]))
        if max(repeat) in translator:
            translator.pop(max(repeat))


def transformMolecules(parser, database, configurationFile, namingConventions,
                       speciesEquivalences=None, bioGridFlag=False):
    """
    main method. Receives a parser configuration, a configurationFile and a
    list of user defined species equivalences and returns a dictionary
    containing an atomized version of the model
    Keywords:
        ---parser: data structure containing the reactions and species we will use
        ---database:data structure containing the result of the outgoing translation
        ---configurationFile
        ---speciesEquivalences:predefined species
    """
    '''
    import cProfile, pstats, StringIO
    pr = cProfile.Profile()
    pr.enable()
    '''

    database = createSpeciesCompositionGraph(parser, database, configurationFile, namingConventions,
                                             speciesEquivalences=speciesEquivalences, bioGridFlag=bioGridFlag)

    for element in database.artificialEquivalenceTranslator:
        if element not in database.eequivalenceTranslator:
            database.eequivalenceTranslator[element] = []
        database.eequivalenceTranslator[element].extend(
            database.artificialEquivalenceTranslator[element])

    # special handling for double modifications like double phosporylation
    # FIXME: this needs to be done in a cleaner way(e.g. getting them
    # from a file instead of being hardcoded)
    doubleModifications = {"Double-Phosporylation": "Phosporylation"}

    for element in doubleModifications:

        if doubleModifications[element] not in database.eequivalenceTranslator:
            continue
        if element not in database.eequivalenceTranslator:
            database.eequivalenceTranslator[element] = []

        baseElements = [
            x[0] for x in database.eequivalenceTranslator[doubleModifications[element]]]
        modifiedElements = [
            x[1] for x in database.eequivalenceTranslator[doubleModifications[element]]]
        #deleteEquivalences = [baseElements.index(x) for x in baseElements if x in modifiedElements]

        deleteEquivalences = [(x, modifiedElements[baseElements.index(x)])
                              for x in baseElements if x in modifiedElements]

        for eq in deleteEquivalences:
            if eq not in database.eequivalenceTranslator[element]:
                database.eequivalenceTranslator[element].append(eq)

        for eq in deleteEquivalences:

            if eq in database.eequivalenceTranslator[doubleModifications[element]]:
                database.eequivalenceTranslator[
                    doubleModifications[element]].remove(eq)

    for modification in database.tmpEquivalence:
        for candidates in database.tmpEquivalence[modification]:
            for instance in candidates:
                addToDependencyGraph(
                    database.eequivalenceTranslator, modification, instance)

    database.weights = sorted(
        database.weights, key=lambda rule: (rule[1], len(rule[0])))

    atomize(database.prunnedDependencyGraph, database.weights, database.translator, database.reactionProperties,
            database.eequivalenceTranslator2, bioGridFlag, database.sbmlAnalyzer, database, parser)

    onlySynDec = len(
        [x for x in database.classifications if x not in ['Generation', 'Decay']]) == 0
    propagateChanges(database.translator, database.prunnedDependencyGraph)

    # sanityCheck(database.translator)
    '''
    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(10)
    print s.getvalue()
    '''
    database.assumptions = assumptions
    return database.translator, onlySynDec
