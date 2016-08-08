# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 17:42:31 2011

@author: proto
"""
from copy import deepcopy, copy
import writer.bnglWriter as writer
log = {'species': [], 'reactions': []}
import re
from collections import Counter
from collections import defaultdict
import math as pymath
from utils.util import logMess, TranslationException
import libsbml

def factorial(x):
    temp = x
    acc = 1
    while temp > 0:
        acc *= temp
        temp -= 1
    return acc


def comb(x, y, exact=True):
    return factorial(x) / (factorial(y) * factorial(x - y))


bioqual = ['BQB_IS', 'BQB_HAS_PART', 'BQB_IS_PART_OF', 'BQB_IS_VERSION_OF',
           'BQB_HAS_VERSION', 'BQB_IS_HOMOLOG_TO',
           'BQB_IS_DESCRIBED_BY', 'BQB_IS_ENCODED_BY', 'BQB_ENCODES', 'BQB_OCCURS_IN',
           'BQB_HAS_PROPERTY', 'BQB_IS_PROPERTY_OF', 'BQB_HAS_TAXON', 'BQB_UNKNOWN']

modqual = ['BQM_IS', 'BQM_IS_DESCRIBED_BY', 'BQM_IS_DERIVED_FROM', 'BQM_IS_INSTANCE_OF', 'BQM_HAS_INSTANCE', 'BQM_UNKNOWN']

annotationHeader = {'BQB':'bqbiol','BQM':'bmbiol'}

def unrollSBMLFunction(function, sbmlFunctions):
    '''
    remove calls to functions inside functions
    '''
    modificationFlag = True
    recursionIndex = 0
    # remove calls to other sbml functions
    while modificationFlag and recursionIndex <20:
        modificationFlag = False
        for sbml in sbmlFunctions:
            if sbml in function:
                temp = writer.extendFunction(function, sbml, sbmlFunctions[sbml])
                if temp != function:
                    function = temp
                    modificationFlag = True
                    recursionIndex +=1
                    break

    return function


class SBML2BNGL:
    '''
    contains methods for extracting and formatting those sbml elements
    that are translatable into bngl
    '''
    def __init__(self, model, useID=True):

        self.useID = useID
        self.model = model
        self.tags = {}
        self.speciesUnits = {}
        self.isConversion = True
        self.boundaryConditionVariables = []
        self.speciesDictionary = {}
        self.speciesMemory = []
        self.speciesAnnotationDict = None
        self.reactionDictionary = {}
        self.speciesAnnotation = None
        self.speciesCompartments = None
        self.unitDefinitions = self.getUnitDefinitions()
        self.convertSubstanceUnits = False

        self.getSpecies()
        
    def setConversion(self,conversion):
        self.isConversion = conversion

    def reset(self):
        self.tags = {}
        self.boundaryConditionVariables = []
        self.speciesDictionary = {}
        self.speciesMemory = []
        self.getSpecies()
        self.reactionDictionary = {}

    def static_var(varname, value):
        def decorate(func):
            setattr(func, varname, value)
            return func
        return decorate

    def extractModelAnnotation(self):
        metaInformation = {}
        annotation = self.model.getAnnotation()
        lista = libsbml.CVTermList()
        libsbml.RDFAnnotationParser.parseRDFAnnotation(annotation, lista)
        for idx in range(lista.getSize()):
            # biol,qual = lista.get(idx).getBiologicalQualifierType(), lista.get(idx).getModelQualifierType()
            qualifierType = lista.get(idx).getQualifierType()
            qualifierDescription = bioqual[lista.get(idx).getBiologicalQualifierType()] if qualifierType \
                else modqual[lista.get(idx).getModelQualifierType()]
            if qualifierDescription not in metaInformation:
                metaInformation[qualifierDescription] = set([])
            for idx2 in range(0, lista.get(idx).getResources().getLength()):
                resource = lista.get(idx).getResources().getValue(idx2)
                metaInformation[qualifierDescription].add(resource)
        return metaInformation

    def getMetaInformation(self, additionalNotes):

        # get unit information
        unitList = self.getUnitDefinitions()

        metaInformation = self.extractModelAnnotation()
        modelHistory = self.model.getModelHistory()
        if modelHistory:
            try:
                tmp = libsbml.ModelHistory.getCreator(self.model.getModelHistory(), 0).getFamilyName()
                tmp += ' ' + libsbml.ModelHistory.getCreator(self.model.getModelHistory(), 0).getGivenName()
                metaInformation['creatorEmail'] = "'" + libsbml.ModelHistory.getCreator(self.model.getModelHistory(), 0).getEmail() + "'"
                metaInformation['creatorName'] = "'" + tmp + "'"
            except:
                metaInformation['creatorEmail'] = "''"
                metaInformation['creatorName'] = "''"

        metaInformation.update(additionalNotes)

        metaString = '###\n'
        for element in metaInformation:
            if type(metaInformation[element]) == set:
                metaInformation[element] = list(metaInformation[element])
            metaString += '#@{0}:{1}\n'.format(element,(metaInformation[element]))
        metaString += '###\n'
        return metaString


    def isSameNameDifferentCompartment(self, name):
        speciesList = []
        for species in self.model.getListOfSpecies():
            if species.getName() == name:
                speciesList.append(species.getCompartment())

        return len(speciesList) == len(set(speciesList))


    def getRawSpecies(self, species,parameters=[], logEntries=True):
        '''
        *species* is the element whose SBML information we will extract
        this method gets information directly
        from an SBML related to a particular species.
        It returns id,initialConcentration,(bool)isconstant and isboundary,
        and the compartment
        It also accounts for the fact that sometimes ppl use the same name for 
        molecules with different identifiers
        '''
        identifier = species.getId()
        name = species.getName()

        if name == '':
            name = identifier
        if species.isSetInitialConcentration():
            initialValue = species.getInitialConcentration()
        else:
            initialValue = species.getInitialAmount()
        isConstant = species.getConstant()
        isBoundary = species.getBoundaryCondition()
        # FIXME: this condition means that a variable/species can be changed
        # by rules and/or events. this means that we effectively need a variable
        # changed by a function that tracks this value, and all references
        # to this observable have to be changed to the referrencing variable.
        # http://sbml.org/Software/libSBML/docs/java-api/org/sbml/libsbml/Species.html
        if isBoundary and not isConstant:
            isConstant = True
            if not species.isSetInitialConcentration() \
                and not species.isSetInitialAmount():
                initialValue = 1
        compartment = species.getCompartment()
        boundaryCondition = species.getBoundaryCondition()
        standardizedName = standardizeName(name)

        # if its a species that appends the compartment name remove it, it is not necessary in bionetgen
        if standardizedName.endswith('_{0}'.format(compartment)):
            standardizedName = standardizedName.replace('_{0}'.format(compartment), '')
        
        #if standardizedName in ['Source','Trash','Sink']:
        #    standardizedName = '0'
        if standardizedName in parameters:
            standardizedName = 'sp_{0}'.format(standardizedName)

        # it cannot start with a number
        if standardizedName[:1].isdigit():
            standardizedName = 's' + standardizedName


        # two species cannot have the same name. Ids are unique but less informative, however typically species can be differentiated
        # by compartment
        if logEntries and standardizedName != '0':

            if standardizedName in self.speciesMemory:
                if len(list(self.model.getListOfCompartments())) == 1:
                    standardizedName += '_' + species.getId()
                else:
                    # we can differentiate by compartment tag, no need to attach it to the name
                    # however an actual check needs to be made to make sure that the modeler is actually 
                    # changing compartment information. If not, use id to differentiate.
                    if not self.isSameNameDifferentCompartment(species.getName()):
                        standardizedName += '_' + species.getId()
            self.speciesMemory.append(standardizedName)

        if boundaryCondition:
            self.boundaryConditionVariables.append(standardizedName)
        self.speciesDictionary[identifier] = standardizedName
        returnID = identifier if self.useID else \
            self.speciesDictionary[identifier]

        if species.isSetSubstanceUnits():
            self.speciesUnits[returnID] = species.getSubstanceUnits()

        values = {}
        values['returnID'] = returnID
        if species.isSetInitialConcentration():
            values['initialConcentration'] = initialValue
            values['initialAmount'] = -1
        elif species.isSetInitialAmount():
            values['initialAmount'] = initialValue
            values['initialConcentration'] = -1
        else:
            values['initialAmount'] = -1
            values['initialConcentration'] = -1

        values['isConstant'] = isConstant
        values['isBoundary'] = isBoundary
        values['compartment'] = compartment
        values['name'] = name
        values['identifier'] = identifier
        return values

    def getPrunnedTree(self, math, remainderPatterns, artificialObservables={}):
        """
        walks through a series of * nodes and removes the remainder reactant factors
        arg:remainderPatterns: argumetns to be removed from the tree
        it also changes references to time variables to the keyword 'Time'

        artificialObservables: species that are changed through an sbml assignment rule. their
        usage in bng requires special handling.
        """
        swapDict = {libsbml.AST_NAME_TIME : 'Time'}
        for node in [x for x in math.getLeftChild(), math.getRightChild() if x != None]:
            if node.getType() in swapDict.keys():
                node.setName(swapDict[node.getType()])
        if math.getCharacter() == '+' and math.getNumChildren() == 1:
            math = math.getLeftChild()
        while (math.getCharacter() == '*' or math.getCharacter() == '/') and len(remainderPatterns) > 0:
            if libsbml.formulaToString(math.getLeftChild()) in remainderPatterns:
                remainderPatterns.remove(libsbml.formulaToString(math.getLeftChild()))
                if math.getCharacter() == '*':
                    math = math.getRightChild()
                else:
                    math.getLeftChild().setValue(1)
            elif libsbml.formulaToString(math.getRightChild()) in remainderPatterns:
                remainderPatterns.remove(libsbml.formulaToString(math.getRightChild()))
                math = math.getLeftChild()
            else:
                if(math.getLeftChild().getCharacter()) == '*':
                    math.replaceChild(0, self.getPrunnedTree(math.getLeftChild(), remainderPatterns))
                if(math.getRightChild().getCharacter()) == '*':
                    math.replaceChild(math.getNumChildren() - 1,self.getPrunnedTree(math.getRightChild(), remainderPatterns))
                break
        return math
        
    def getIsTreeNegative(self, math):
        '''
        walks through a series of * nodes and detects whether there's a negative factor
        fixme: we should actually test if the number of negative factors is odd
        right now we are relying on  the modelers not being malicious
        when writing their rates laws.
        '''

        if (math.getCharacter() == '*' or math.getCharacter() == '/'):
            if math.getLeftChild().isUMinus():
                math.setCharacter('+')
                #math.getLeftChild().setValue(long(libsbml.formulaToString(math.getLeftChild())[1:]))
                return True
            elif math.getLeftChild().getNumChildren() == 0 and libsbml.formulaToString(math.getLeftChild()).startswith('-'):
                math.getLeftChild().setValue(long(libsbml.formulaToString(math.getLeftChild())[1:]))
                return True
            elif math.getRightChild().isUMinus():
                math.setCharacter('+')
                #math.getRightChild().setValue(long(libsbml.formulaToString(math.getRightChild())[1:]))
                return True
            elif math.getRightChild().getNumChildren() == 0 and libsbml.formulaToString(math.getRightChild()).startswith('-'):
                math.getRightChild().setValue(long(libsbml.formulaToString(math.getRightChild())[1:]))
                return True

            else:
                if(math.getLeftChild().getCharacter()) in ['*', '/', '-']:
                    if self.getIsTreeNegative(math.getLeftChild()):
                        return True
                if(math.getRightChild().getCharacter()) in ['*', '/', '-']:
                    if self.getIsTreeNegative(math.getRightChild()):
                        return True
        elif math.getCharacter() == '-' and math.getNumChildren() == 1:
            math.setCharacter('+')
            return True
        return False

    def getUnitDefinitions(self):
        self.unitDictionary = {}
        for unitDefinition in self.model.getListOfUnitDefinitions():
            unitList = []
            for unit in unitDefinition.getListOfUnits():
                correctedUnit = libsbml.Unit.convertToSI(unit)
                #unitList.append({'kind':unit.getKind(), 'scale':unit.getScale(),'multiplier':unit.getMultiplier(), 'exponent': unit.getExponent(), 'name':name})
                for unit2 in correctedUnit.getListOfUnits():
                    name = unit.getName() if unit.getName() else unitDefinition.getName()
                    unitList.append({'kind':unit2.getKind(), 'scale': unit2.getScale(), 'multiplier':unit2.getMultiplier(), 'exponent': unit2.getExponent(), 'name':name})
            self.unitDictionary[unitDefinition.getId()] = unitList
        return self.unitDictionary

    def preProcessStoichiometry(self, reactants):
        '''
        checks for reactants with the same name in the reactant list. This
        is mainly to account for reactants that dont have the stoichiometry
        flag properly set and instead appear repeated
        '''
        uniqueReactantDict = defaultdict(int)
        for reactant in reactants:
            uniqueReactantDict[reactant[0]] += reactant[1]
        return [(x, uniqueReactantDict[x]) for x in uniqueReactantDict]

    def removeFactorFromMath(self, math, reactants, products, artificialObservables):
        '''
        it also adds symmetry factors. this checks for symmetry in the species names
        s

        artificialObservables: species names that are changed through assignment rules. their use requires special care when calculating a rate
        '''
        ifStack = Counter()
        remainderPatterns = []
        highStoichoiMetryFactor = 1
        processedReactants = self.preProcessStoichiometry(reactants)
        for x in processedReactants:
            highStoichoiMetryFactor *= factorial(x[1])
            y = [i[1] for i in products if i[0] == x[0]]
            y = y[0] if len(y) > 0 else 0
            # TODO: check if this actually keeps the correct dynamics
            # this is basically there to address the case where theres more products
            # than reactants (synthesis)
            if x[1] > y:
                highStoichoiMetryFactor /= comb(int(x[1]), int(y), exact=True)
            for counter in range(0, int(x[1])):
                remainderPatterns.append(x[0])
        
        #for x in products:
        #    highStoichoiMetryFactor /= math.factorial(x[1])
        #remainderPatterns = [x[0] for x in reactants]
        #print remainderPatterns, artificialObservables.keys()
        math = self.getPrunnedTree(math, remainderPatterns)

        rateR = libsbml.formulaToString(math)

        for element in remainderPatterns:
            ifStack.update([element])
        for element in ifStack:
            if ifStack[element] > 1:
                rateR = 'if({0}>0, {1}/({0}^{2}),0)'.format(element, rateR, ifStack[element])
            else:
                rateR = 'if({0}>0, {1}/{0},0)'.format(element, rateR)

        numFactors = max(math.getNumChildren(), len(ifStack))
        if pymath.isinf(highStoichoiMetryFactor):
            rateR = '{0} * 1e20'.format(rateR)
            logMess('ERROR:SIM204','Found usage of "inf" inside function {0}'.format(rateR))
        elif highStoichoiMetryFactor != 1:
            rateR = '{0} * {1}'.format(rateR, int(highStoichoiMetryFactor))
            # we are adding a factor to the rate so we need to account for it when 
            # we are constructing the bngl equation (we dont want constrant expressions in there)
            numFactors = max(1, numFactors)
        return rateR, numFactors

    def isAmount(self, reactantName):
        for species in self.model.getListOfSpecies():
            if species.getName() == reactantName:
                if species.isSetInitialAmount():
                    return True
        return False


    def analyzeReactionRate(self, math, compartmentList, reversible, rReactant, rProduct, reactionID, parameterFunctions, rModifier=[], sbmlFunctions={}):
        """
        This functions attempts to obtain the left and right hand sides of a rate reaction 
        function given a MathML tree. It also removes compartments and chemical factors from the function

        Keyword arguments:
            math -- the MathML math object
            compartmentList -- a list of all the compartments in the system
            reversible -- a boolean indicated whether there's should be a backward rate
            rReactant -- a string list of the reactants.
            rProduct -- a string list of the products.
            sbmlFunctions -- a list of possible nested functiosn taht we need to remove
        """
        removedCompartments = copy(compartmentList)
        math = self.getPrunnedTree(math, compartmentList)

        mathstring = libsbml.formulaToString(math)
        #unroll sbml functions
        mathstring =  unrollSBMLFunction(mathstring, sbmlFunctions)

        temp = libsbml.parseFormula(mathstring)
        # if we could parse back after all those modifications...
        if temp:
            math = temp

        moleFlag = False
        if any([element[0] in self.speciesUnits for element in rReactant]):
            for element in rReactant:
                if element[0] in self.speciesUnits and self.speciesUnits[element[0]] in self.unitDefinitions:
                    if self.unitDefinitions[self.speciesUnits[element[0]]] == 23:
                        moleFlag = True
                    else:
                        moleFlag = False
        else:
            if 'substance' not in self.unitDefinitions:
                moleFlag = True
            elif any([x['kind'] == 23 for x in self.unitDefinitions['substance']]):
                moleFlag = True
        #else if all([x['name'] != 'units' for x in self.unitDefinitions]):
        #    print self.unitDefinitions['substance']
        #divide by avogadros number to get volume per number per second units

        removedCompartments = [x for x in removedCompartments if x not in compartmentList]
        if reversible:
            if math.getCharacter() == '-':
                if math.getNumChildren() > 1:
                    rateL, nl = (self.removeFactorFromMath(
                                 math.getLeftChild().deepCopy(), rReactant,
                                 rProduct, parameterFunctions))
                    rateR, nr = (self.removeFactorFromMath(
                                 math.getRightChild().deepCopy(), rProduct,
                                 rReactant, parameterFunctions))
                else:
                    rateR, rateL, nr, nl = self.analyzeReactionRate(math.getChild(0), compartmentList,
                                                                    reversible, rProduct, rReactant, reactionID, parameterFunctions, rModifier, sbmlFunctions)
            elif math.getCharacter() == '+' and math.getNumChildren() > 1:
                if(self.getIsTreeNegative(math.getRightChild())):
                    rateL, nl = (self.removeFactorFromMath(
                    math.getLeftChild().deepCopy(), rReactant, rProduct, parameterFunctions))
                    rateR, nr = (self.removeFactorFromMath(
                    math.getRightChild().deepCopy(), rProduct, rReactant, parameterFunctions))
                elif(self.getIsTreeNegative(math.getLeftChild())):
                    rateR, nr = (self.removeFactorFromMath(
                                 math.getLeftChild().deepCopy(), rProduct, rReactant, parameterFunctions))
                    rateL, nl = (self.removeFactorFromMath(
                    math.getRightChild().deepCopy(), rReactant, rProduct, parameterFunctions))
                else:
                    rateL, nl = self.removeFactorFromMath(math.deepCopy(), rReactant,
                                                          rProduct, parameterFunctions)
                    rateL = "if({0}>= 0,{0},0)".format(rateL)
                    rateR, nr = self.removeFactorFromMath(math.deepCopy(), rProduct,
                                                          rReactant, parameterFunctions)
                    rateR = "if({0}< 0,-({0}),0)".format(rateR)
                    nl, nr = 1, 1

            else:
                # reaction is bidirectional but i can't separate function into
                # left hand side and right hand side
                rateL, nl = self.removeFactorFromMath(math.deepCopy(), rReactant,
                                                      rProduct, parameterFunctions)
                rateR, nr = self.removeFactorFromMath(math.deepCopy(), rProduct,
                                                      rReactant, parameterFunctions)
                if nl > 0:
                    if nr == 0 and rateR not in parameterFunctions:
                        rateL = '0'
                        nl = -1
                        logMess('INFO:SIM001', 'In reaction {0}, the left hand side has been determined \
to never activate and has been to rate 0'.format(reactionID))
                    else:
                        rateL = "if({0} >= 0, {0}, 0)".format(rateL)
                        nl = 1
                if nr > 0:
                    if nl == 0 and rateL not in parameterFunctions:
                        rateR = '0'
                        nr = -1
                        logMess('INFO:SIM002', 'In reaction {0}, the right hand side has been determined \
to never activate (rate is never negative), setting reaction to unidirectional'.format(reactionID))
                    else:
                        rateR = "if({0} < 0, -({0}), 0)".format(rateR)
                        nr = 1
                if ((nl == 0 and nr > 0) or (nr == 0 and nl > 0)) and (rateL in parameterFunctions or rateR in parameterFunctions):
                    logMess('WARNING:SIM102', 'In reaction {0}, rates cannot be divided into left hand side and right hand side \
but reaction is marked as reversible'.format(reactionID))

                #nl, nr = 1,1
        else:
            rateL, nl = (self.removeFactorFromMath(math.deepCopy(),
                                                   rReactant, rProduct, parameterFunctions))

            rateR, nr = '0', '-1'


            
        #cBNGL and SBML treat the behavior of compartments in rate laws differently so we have to compensate for that
        if len(removedCompartments) > 0:
            #if the species initial conditions were defined as concentrations then correct for it and transform it to absolute counts
            if len(rReactant) == 2 and not (self.isAmount(rReactant[0][0]) or self.isAmount(rReactant[1][0])):
                if moleFlag:
                    rateL = '({0}) / 6.022e23'.format(rateL)
                    nl += 1
                #rateL = '({0}) * ({1})'.format(rateL,' * '.join(removedCompartments))
                #nl += 1
                #pass
            elif len(rModifier) > 0:
                pass
                #rateL = '({0}) / ({1})'.format(rateL,' * '.join(removedCompartments))
                #nl += 1
                #pass
            if nr != '-1':
                if len(rProduct) == 2 and len(rReactant) == 1 and not (self.isAmount(rProduct[0][0]) or self.isAmount(rProduct[1][0])):
                    if moleFlag:
                        rateR = '({0}) / 6.022e23'.format(rateR)
                        nl += 1

                    #rateR = '({0}) * {1}'.format(rateR,' * '.join(removedCompartments))
                    #nr += 1
                    #pass
        return rateL, rateR, nl, nr

    def __getRawRules(self, reaction, symmetryFactors, parameterFunctions, translator, sbmlfunctions):

        zerospecies = ['emptyset','trash','sink','source']
        if self.useID:
            reactant = [(reactant.getSpecies(), reactant.getStoichiometry())
                        for reactant in reaction.getListOfReactants() if
                        reactant.getSpecies().lower() not in zerospecies and reactant.getStoichiometry() not in [0,'0']]
            product = [(product.getSpecies(), product.getStoichiometry())
                       for product in reaction.getListOfProducts() if product.getSpecies().lower()
                       not in zerospecies and product.getStoichiometry() not in [0,'0']]
        else:
            reactant = [(self.speciesDictionary[rElement.getSpecies()], rElement.getStoichiometry(), rElement.getSpecies())
                        for rElement in reaction.getListOfReactants() if self.speciesDictionary[rElement.getSpecies()].lower() not in zerospecies
                        and rElement.getStoichiometry() not in [0,'0']]
            product = [(self.speciesDictionary[rProduct.getSpecies()], rProduct.getStoichiometry(), rProduct.getSpecies())
                       for rProduct in reaction.getListOfProducts() if self.speciesDictionary[rProduct.getSpecies()].lower() not in zerospecies
                       and rProduct.getStoichiometry() not in [0,'0']]
        kineticLaw = reaction.getKineticLaw()
        reversible = reaction.getReversible()
        if kineticLaw is None:
            return {'reactants': reactant, 'products': product, 'parameters': [], 'rates': ['0', '0'],
                    'reversible': reversible, 'reactionID': reaction.getId(), 'numbers': [0, 0], 'modifiers': None}

        rReactant = []
        rProduct = []

        # in case a given species was defined as the zero molecule don't include it in the rate correction method
        for x in reaction.getListOfReactants():
            if x.getSpecies().lower() not in zerospecies and x.getStoichiometry() not in [0, '0']:
                if not x.getConstant() and pymath.isnan(x.getStoichiometry()):
                    logMess("ERROR:SIM241", "BioNetGen does not support non constant stoichiometries. Reaction {0} is not correctly translated".format(reaction.getId()))
                    raise TranslationException(reaction.getId())
                else:
                    speciesName = self.speciesDictionary[x.getSpecies()] if x.getSpecies() in self.speciesDictionary else ''
                    if speciesName in translator and str(translator[speciesName]) == '0':
                        continue
                    rReactant.append((x.getSpecies(), x.getStoichiometry()))
        for x in reaction.getListOfProducts():
            if x.getSpecies().lower() not in zerospecies and x.getStoichiometry() not in [0, '0']:
                if not x.getConstant() and pymath.isnan(x.getStoichiometry()):
                    logMess("ERROR:SIM241", "BioNetGen does not support non constant stoichiometries. Reaction {0} is not correctly translated".format(reaction.getId()))
                    raise TranslationException(reaction.getId())
                else:
                    speciesName = self.speciesDictionary[x.getSpecies()] if x.getSpecies() in self.speciesDictionary else ''
                    if speciesName in translator and str(translator[speciesName]) == '0':
                        continue
                    rProduct.append((x.getSpecies(), x.getStoichiometry()))
        #rReactant = [(x.getSpecies(), x.getStoichiometry()) for x in reaction.getListOfReactants() if x.getSpecies() not in ['EmptySet']]
        #rProduct = [(x.getSpecies(), x.getStoichiometry()) for x in reaction.getListOfProducts() if x.getSpecies() not in ['EmptySet']]
        rModifiers = [x.getSpecies() for x in reaction.getListOfModifiers() if x.getSpecies() != 'EmptySet']
        parameters = [(parameter.getId(), parameter.getValue(), parameter.getUnits()) for parameter in kineticLaw.getListOfParameters()]

        rateL = rateR = nl = nr = None
        if True:
            # TODO: For some reason creating a deepcopy of this screws everything up, even
            # though its what we should be doing
            # update: apparently the solution was to use copy instead of deepcopy. This is because
            # the underlying swig code in c was causing conflicts when copied. make sure this actually works
            math = copy(kineticLaw.getMath())
            math = math.deepCopy()
            # get a list of compartments so that we can remove them
            compartmentList = []
            for compartment in (self.model.getListOfCompartments()):
                compartmentList.append(compartment.getId())

            # remove compartments from expression. also separate left hand and right hand side


            rateL, rateR, nl, nr = self.analyzeReactionRate(math, compartmentList,
                reversible, rReactant, rProduct, reaction.getId(), parameterFunctions, rModifiers, sbmlfunctions)

            if rateR == '0':
                reversible = False
            if symmetryFactors[0] > 1:
                rateL = '({0})/{1}'.format(rateL, symmetryFactors[0])
            if symmetryFactors[1] > 1:
                rateR = '({0})/{1}'.format(rateR, symmetryFactors[1])
            if not self.useID:
                rateL = self.convertToName(rateL)
                rateR = self.convertToName(rateR)
            if reversible:
                pass
            # return compartments if the reaction is unimolecular
            # they were removed in the first palce because its easier to handle
            # around the equation in tree form when it has less terms
            '''
            if len(self.model.getListOfCompartments()) > 0:
                for compartment in (self.model.getListOfCompartments()):
                    if compartment.getId() not in compartmentList:
                        if len(rReactant) != 2:
                            rateL = '{0} * {1}'.format(rateL,compartment.getSize())
                        if len(rProduct) != 2:
                             rateR = '{0} * {1}'.format(rateR,compartment.getSize())
            '''
        return {'reactants': reactant, 'products': product, 'parameters': parameters, 'rates': [rateL, rateR],
                'reversible': reversible, 'reactionID': reaction.getId(), 'numbers': [nl, nr], 'modifiers': rModifiers}

    def getReactionCenter(self, reactant, product, translator):
        rcomponent = Counter()
        pcomponent = Counter()

        for element in reactant:
            if element[0] in translator:
                for molecule in translator[element[0]].molecules:
                    for component in molecule.components:
                        molecule.sort()
                        rcomponent.update(Counter([(molecule.name, component.name, len(component.bonds) > 0, component.activeState)]))
        for element in product:
            if element[0] in translator:
                for molecule in translator[element[0]].molecules:
                    molecule.sort()
                    for component in molecule.components:
                        molecule.sort()
                        pcomponent.update(Counter([(molecule.name, component.name, len(component.bonds) > 0, component.activeState)]))
        reactionCenter = [(x[0], x[1]) for x in rcomponent for y in pcomponent if (x[0], x[1]) == (y[0], y[1]) and x != y and rcomponent[x] != pcomponent[y]]
        rreactionCenter = [(x[0], x[1]) for x in pcomponent for y in rcomponent if (x[0], x[1]) == (y[0], y[1]) and x != y and pcomponent[x] != rcomponent[y]]
        return reactionCenter, rreactionCenter

    def updateComponentCount(self, counterArray, reference, updateValue):
        for element in counterArray:
            if reference in counterArray[element]:
                counterArray[element][reference] += updateValue

    def reduceComponentSymmetryFactors(self, reaction, translator, functions):
        '''
        create symmetry factors for reactions with components and species with
        identical names. This checks for symmetry in the components names then.
        '''
        zerospecies = ['emptyset','trash','sink','source']
        if self.useID:
            reactant = [(rElement.getSpecies(), rElement.getStoichiometry())
            for rElement in reaction.getListOfReactants() if
            rElement.getSpecies() != 'EmptySet']
            product = [(product.getSpecies(), product.getStoichiometry())
            for product in reaction.getListOfProducts() if product.getSpecies()
            != 'EmptySet']
        else:
            reactant = [(self.speciesDictionary[rElement.getSpecies()], rElement.getStoichiometry()) for rElement in reaction.getListOfReactants()]
            product = [(self.speciesDictionary[rProduct.getSpecies()], rProduct.getStoichiometry()) for rProduct in reaction.getListOfProducts()]
        kineticLaw = reaction.getKineticLaw()
        reversible = reaction.getReversible()

        if kineticLaw is None:
            return 1, 1
        rReactant = rProduct = []
        
        for x in reaction.getListOfReactants():
            if x.getSpecies().lower() not in zerospecies \
                        and x.getStoichiometry() not in [0, '0'] and pymath.isnan(x.getStoichiometry()):
                if not x.getConstant():
                    logMess("ERROR:SIM241", "BioNetGen does not support non constant stoichiometries. Reaction {0} is not correctly translated".format(reaction.getId()))
                    return 1, 1
                else:
                    rReactant.append(x.getSpecies(), x.getStoichiometry())

        for x in reaction.getListOfProducts():
            if x.getSpecies().lower() not in zerospecies \
                        and x.getStoichiometry() not in [0, '0'] and pymath.isnan(x.getStoichiometry()):
                if not x.getConstant():
                    logMess("ERROR:SIM241", "BioNetGen does not support non constant stoichiometries. Reaction {0} is not correctly translated".format(reaction.getId()))
                    return 1, 1
                else:
                    rProduct.append(x.getSpecies(), x.getStoichiometry())
        
        
        # TODO: For some reason creating a deepcopy of this screws everything up, even
        # though its what we should be doing
        rcomponent = defaultdict(Counter)
        pcomponent = defaultdict(Counter)
        
        #get the total count of components in the reactants and products
        #e.g. components across diffent species
        freactionCenter,breactionCenter = self.getReactionCenter(reactant, product, translator)
        
        for element in reactant:
            if element[0] in translator:
                
                for molecule in translator[element[0]].molecules:
                    for component in molecule.components:
                        molecule.sort()
                        componentList = Counter([(molecule.signature(freactionCenter))])
                        for _ in range(0,int(element[1])):
                            rcomponent[(molecule.name,component.name,len(component.bonds)>0,component.activeState)].update(componentList)
                        

        
        for element in product:
            if element[0] in translator:
                for molecule in translator[element[0]].molecules:
                    molecule.sort()
                    for component in molecule.components:

                        componentList = Counter([(molecule.signature(breactionCenter))])
                        for _ in range(0,int(element[1])):
                            pcomponent[(molecule.name,component.name,len(component.bonds)>0,component.activeState)].update(componentList)

        '''
        only keep information for reaction centers
        '''
        reactionCenters = [(x[0],x[1]) for x in rcomponent for y in pcomponent if (x[0], x[1]) == (y[0], y[1]) and x != y]    
        rcomponent= {x:rcomponent[x] for x in rcomponent if (x[0], x[1]) in reactionCenters}
        pcomponent= {x:pcomponent[x] for x in pcomponent if (x[0], x[1]) in reactionCenters}

        #is the number of components across products and reactants the same?
        #eg is there any DeleteMolecules action
        pcorrectionFactor = 1
        rcorrectionFactor = 1
        rStack = []
        pStack = []
        '''
        if a reaction can take place in several ways account for it in the reaction 
        rate (this is specially important in dimer and trimer binding)
        pcomponent[element] < rcomponent[element] asks if an specific instance
        of a component decreases in number from a reactant to a product
        for example if there are 3 A(b)'s and one binds, we will have 2 A(b)'s
        in the product  
        '''
        rcomponentTemp = deepcopy(rcomponent)
        pcomponentTemp = deepcopy(pcomponent)
        
        #calculate actual symmetry factors
        for key in rcomponent:
            if key in pcomponent:
                for element in rcomponent[key]:
                    if rcomponent[key] ==1:
                        continue
                    #if theres a component on one side of the equation that
                    #appears a different number of times on the other side of the equation
                    if element in pcomponent[key]:
                        if pcomponent[key][element] < rcomponent[key][element] and set([key[0].lower(),key[1].lower()]) not in rStack:
                            rcorrectionFactor *= comb(rcomponent[key][element],pcomponent[key][element],exact=1)
                            rStack.append(set([key[0].lower(),key[1].lower()]))
                            #once we choose a component for a previous action
                            #this limits the options for subsequent actions
                            #although substracting one from the target sites
                            #may not be the right option. double check.
                            self.updateComponentCount(pcomponent,element,-1)
                    else:
                        for element2 in pcomponent[key]:
                            if pcomponent[key][element2] < rcomponent[key][element] and set([key[0].lower(),key[1].lower()]) not in rStack:
                                rcorrectionFactor *= comb(rcomponent[key][element],pcomponent[key][element2],exact=1)
                                rStack.append(set([key[0].lower(),key[1].lower()]))
                                self.updateComponentCount(pcomponent,element2,-1)

        rcomponent = rcomponentTemp
        pcomponent = pcomponentTemp
        
        if reversible:
            for key in pcomponent:
                if key in rcomponent:
                    for element in pcomponent[key]:
                        if pcomponent[key] ==1:
                            continue
                        if element in rcomponent[key]:
                            if rcomponent[key][element] < pcomponent[key][element] and set([key[0].lower(),key[1].lower()]) not in pStack:
                                pcorrectionFactor *= comb(pcomponent[key][element],rcomponent[key][element],exact=1)
                                pStack.append(set([key[0].lower(),key[1].lower()]))
                                self.updateComponentCount(rcomponent,element,-1)

                        else:
                            for element2 in rcomponent[key]:
                                if rcomponent[key][element2] < pcomponent[key][element] and set([key[0].lower(),key[1].lower()]) not in pStack:
                                    pcorrectionFactor *= comb(pcomponent[key][element],rcomponent[key][element2],exact=1)
                                    pStack.append(set([key[0].lower(),key[1].lower()]))
                                    self.updateComponentCount(rcomponent,element2,-1)


        return rcorrectionFactor,pcorrectionFactor
        
    def convertToName(self, rate):
        for element in sorted(self.speciesDictionary, key=len, reverse=True):
            while True:
                oldRate = rate
                if element in rate:
                    rate = re.sub(r'(\W|^)({0})(\W|$)'.format(element),
                                  r'\1{0}\3'.format(
                                  self.speciesDictionary[element]), rate)
                    if rate != oldRate:
                        continue
                break
            #rate = rate.replace(element,self.speciesDictionary[element])
        return rate

    def __getRawCompartments(self, compartment):
        '''
        Private method used by the getCompartments method 
        '''
        idid = compartment.getId()
        name = compartment.getName()
        size = compartment.getSize()
        dimensions = compartment.getSpatialDimensions()
        if dimensions in [0, 1]:
            logMess('WARNING:SIM103', '{1}-D compartments are not supported. Changing for 2-D compartments for {0}. Please verify this does not affect simulation'.format(name, dimensions))
            dimensions = 2
        #if size != 1:
        #    print '!',
        #return name,3,size
        return idid, dimensions, size, name
        
    def __getRawFunctions(self,function):
        math= function[1].getMath()
        name = function[1].getId()
        
        return name,libsbml.formulaToString(math)

    def getSBMLFunctions(self):
        functions = {}
        for function in enumerate(self.model.getListOfFunctionDefinitions()):
            functionInfo = self.__getRawFunctions(function)
            functions[functionInfo[0]] = (writer.bnglFunction(functionInfo[1],functionInfo[0],[],reactionDict=self.reactionDictionary))
        return functions
            
    def getCompartments(self):
        '''
        Returns an array of triples, where each triple is defined as
        (compartmentName,dimensions,size)
        '''
        compartments = []
        unitDefinitions = self.getUnitDefinitions()
        if 'volume' in unitDefinitions:
            compartments.append('#volume units: {0}'.format('*'.join([x['name'] for x in unitDefinitions['volume']])))
        else:
            compartments.append('#volume units: L')
        for _,compartment in enumerate(self.model.getListOfCompartments()):
            compartmentInfo = self.__getRawCompartments(compartment)
            name = 'cell' if compartmentInfo[0] == '' else compartmentInfo[0]
            if name != compartmentInfo[3]:
                compartments.append("%s  %d  %s #%s" % (name, compartmentInfo[1], compartmentInfo[2], compartmentInfo[3]))
            else:
                compartments.append("%s  %d  %s" % (name, compartmentInfo[1], compartmentInfo[2]))
        return compartments



    def updateFunctionReference(self, reaction, updatedReferences):
        newRate = reaction[3]
        for reference in updatedReferences:
            newRate = re.sub(r'(\W|^)({0})(\W|$)'.format(reference), r'\1{0}\3'.format(updatedReferences[reference]), newRate)

        return newRate

    def getReactions(self, translator={}, isCompartments=False, extraParameters={}, atomize=False, parameterFunctions={}, database= None):
        '''
        @returns: a triple containing the parameters,reactions,functions
        '''

        # @FIXME:this part of the code is there so that we only generate the functions list once through different
        # iterations of this call. This is because we cannot create a clone of the 'math' object for this
        # reaction and it is being permanently changed every call. It's ugly but it works. Change for something
        # better when we figure out how to clone the math object
        if not hasattr(self.getReactions, 'functionFlag'):
            self.getReactions.__func__.functionFlag = False or (not atomize)

        reactions = []
        reactionStructure = []
        parameters = []
        functions = []
        functionTitle = 'functionRate'

        self.unitDefinitions = self.getUnitDefinitions()
        database.rawreactions = []
        if len(self.model.getListOfReactions()) == 0:
            logMess('WARNING:SIM104', 'Model contains no natural reactions, all reactions are produced by SBML rules')
        for index, reaction in enumerate(self.model.getListOfReactions()):
            parameterDict = {}
            # symmetry factors for components with the same name
            sl, sr = self.reduceComponentSymmetryFactors(reaction, translator, functions)
            sbmlfunctions = self.getSBMLFunctions()

            try:
                rawRules = self.__getRawRules(reaction, [sl, sr], parameterFunctions, translator, sbmlfunctions)
                database.rawreactions.append(rawRules)
            except TranslationException as e:
                if(database and database.ignore):
                    reactions.append('#Reaction {0} is not correctly translated. check log for details'.format(e.value))
                    continue
                else:
                    raise TranslationException(e.value + " during reaction processing")

            if len(rawRules['parameters']) > 0:
                for parameter in rawRules['parameters']:
                    """
                    if self.convertSubstanceUnits:
                        #newParameter = self.convertToStandardUnits(parameter[1],self.unitDictionary)
                        parameters.append('r%d_%s %f' % (index+1, parameter[0], parameter[1]))
                        parameterDict[parameter[0]] = parameter[1]
                    else:
                    """
                    parameters.append('r%d_%s %f' % (index + 1, parameter[0], parameter[1]))
                    parameterDict[parameter[0]] = parameter[1]
            compartmentList = [['cell', 1]]
            compartmentList.extend([[self.__getRawCompartments(x)[0], self.__getRawCompartments(x)[2]] for x in self.model.getListOfCompartments()])
            threshold = 0
            if rawRules['numbers'][0] > threshold  or rawRules['rates'][0] in translator:
                functionName = '%s%d()' % (functionTitle, index)
            else:
                # append reactionNumbers to parameterNames
                finalString = str(rawRules['rates'][0])
                for parameter in parameterDict:
                    finalString = re.sub(r'(\W|^)({0})(\W|$)'.format(parameter),
                                         r'\1{0}\3'.format('r{0}_{1}'.format(index + 1, parameter)),
                                         finalString)
                functionName = finalString
            if self.getReactions.functionFlag and 'delay' in rawRules['rates'][0]:
                logMess('ERROR:SIM202', 'BNG cannot handle delay functions in function %s' % functionName)
            if rawRules['reversible']:
                if rawRules['numbers'][0] > threshold or rawRules['rates'][0] in translator:
                    if self.getReactions.functionFlag:
                        functions.append(writer.bnglFunction(rawRules['rates'][0], functionName, rawRules['reactants'], compartmentList, parameterDict, self.reactionDictionary))
                if rawRules['numbers'][1] > threshold  or rawRules['rates'][1] in translator:
                    functionName2 = '%s%dm()' % (functionTitle, index)
                    if self.getReactions.functionFlag:
                        functions.append(writer.bnglFunction(rawRules['rates'][1], functionName2, rawRules['products'],
                                         compartmentList, parameterDict, self.reactionDictionary))
                    self.reactionDictionary[rawRules['reactionID']] = '({0} - {1})'.format(functionName, functionName2)
                    functionName = '{0},{1}'.format(functionName, functionName2)
                else:
                    finalString = str(rawRules['rates'][1])
                    for parameter in parameterDict:
                        finalString = re.sub(r'(\W|^)({0})(\W|$)'.format(parameter),
                                             r'\1{0}\3'.format('r{0}_{1}'.format(index + 1, parameter)),
                                             finalString)
                    functionName = '{0},{1}'.format(functionName, finalString)

            else:

                if rawRules['numbers'][0] > threshold or rawRules['rates'][0] in translator:

                    if self.getReactions.functionFlag:
                        functions.append(writer.bnglFunction(rawRules['rates'][0], functionName, rawRules['reactants'],
                                                             compartmentList, parameterDict, self.reactionDictionary))
                    self.reactionDictionary[rawRules['reactionID']] = '{0}'.format(functionName)
            #reactants = [x for x in rawRules[0] if x[0] not in self.boundaryConditionVariables]
            #products = [x for x in rawRules[1] if x[0] not in self.boundaryConditionVariables]
            reactants = [x for x in rawRules['reactants']]
            products = [x for x in rawRules['products']]
            modifierComment = '#Modifiers({0})'.format(', '.join(rawRules['modifiers'])) if rawRules['modifiers'] else ''
            reactions.append(writer.bnglReaction(reactants, products, functionName, self.tags, translator,
                             (isCompartments or ((len(reactants) == 0 or len(products) == 0) and self.getReactions.__func__.functionFlag)),
                             rawRules['reversible'], reactionName=rawRules['reactionID'], comment=modifierComment))

        if atomize:
            self.getReactions.__func__.functionFlag = True
        return parameters, reactions, functions

    def __getRawAssignmentRules(self, arule):
        variable = arule.getVariable()
        
        # try to separate into positive and negative sections
        if arule.getMath().getCharacter() == '-' and arule.getMath().getNumChildren() > 1 and not arule.isAssignment():
            rateL = libsbml.formulaToString(arule.getMath().getLeftChild())
            if(arule.getMath().getRightChild().getCharacter()) == '*':
                if libsbml.formulaToString(arule.getMath().getRightChild().getLeftChild()) == variable:
                    rateR = libsbml.formulaToString(arule.getMath().getRightChild().getRightChild())
                elif libsbml.formulaToString(arule.getMath().getRightChild().getRightChild()) == variable:
                    rateR = libsbml.formulaToString(arule.getMath().getRightChild().getLeftChild())
                else:
                    rateR = 'if({0}>0,({1})/{0},0)'.format(variable, libsbml.formulaToString(arule.getMath().getRightChild()))
            else:
                rateR = 'if({0}>0,({1})/{0},0)'.format(variable, libsbml.formulaToString((arule.getMath().getRightChild())))
        else:
            rateL = libsbml.formulaToString(arule.getMath())
            rateR = '0'
        if not self.useID:
            rateL = self.convertToName(rateL)
            rateR = self.convertToName(rateR)
            #variable = self.convertToName(variable).strip()
        #print arule.isAssignment(),arule.isRate()
        return variable,[rateL, rateR], arule.isAssignment(), arule.isRate()
        
    def getAssignmentRules(self, zparams, parameters, molecules, observablesDict, translator):
        '''
        this method obtains an SBML rate rules and assignment rules. They
        require special handling since rules are often both defined as rules 
        and parameters initialized as 0, so they need to be removed from the parameters list
        '''

        compartmentList = [['cell',1]]
        compartmentList.extend([[self.__getRawCompartments(x)[0], self.__getRawCompartments(x)[2]] for x in self.model.getListOfCompartments()])

        arules = []
        aParameters = {}
        zRules = zparams
        removeParameters = []
        artificialReactions = []
        artificialObservables = {}
        nonamecounter = 0
        for arule in self.model.getListOfRules():
            rawArule = self.__getRawAssignmentRules(arule)

            #rule has no name
            if rawArule[0] == '':
                logMess('ERROR:SIM215','atomizer has found an sbml rule without a name. {0}'.format(rawArule[1:]))
                rawArule = list(rawArule)
                rawArule[0] = 'noname{0}'.format(nonamecounter)
                nonamecounter += 1
            #tmp.remove(rawArule[0])
            #newRule = rawArule[1].replace('+',',').strip()
            if rawArule[3] == True:
                #it is a rate rule

                if rawArule[0] in self.boundaryConditionVariables:
                    logMess('WARNING:SIM105','rate rules ({0}) \
                    are not properly supported in BioNetGen simulator'.format(rawArule[0]))

                    #aParameters[rawArule[0]] = 'arj' + rawArule[0] 
                    #tmp = list(rawArule)
                    #tmp[0] = 'arj' + rawArule[0]
                    #rawArule = tmp


                rateLaw1 = rawArule[1][0]
                rateLaw2 = rawArule[1][1]
                arules.append(writer.bnglFunction(rateLaw1, 'arRate{0}'.format(rawArule[0]),[],compartments=compartmentList, reactionDict=self.reactionDictionary))
                arules.append(writer.bnglFunction(rateLaw2, 'armRate{0}'.format(rawArule[0]),[],compartments=compartmentList, reactionDict=self.reactionDictionary))
                #moleculeString = str(translator[rawArule[0]]) if rawArule[0] in translator else rawArule[0]
                artificialReactions.append(writer.bnglReaction([], [[self.convertToName(rawArule[0]).strip(),1, rawArule[0]]],'{0},{1}'.format('arRate{0}'.format(rawArule[0]), 'armRate{0}'.format(rawArule[0])), self.tags, translator, isCompartments=True, comment = '#rateLaw'))
                #arules.append(writer.bnglFunction('({0}) - ({1})'.format(rawArule[1][0],rawArule[1][1]), '{0}'.format(rawArule[0]),[],compartments=compartmentList, reactionDict=self.reactionDictionary))
                if rawArule[0] in zparams:
                    removeParameters.append('{0} 0'.format(rawArule[0]))
                    zRules.remove(rawArule[0])
                else:
                    for element in parameters:
                        #TODO: if for whatever reason a rate rule
                        #was defined as a parameter that is not 0
                        #remove it. This might not be exact behavior
                        if re.search('^{0}\s'.format(rawArule[0]), element):
                            logMess("WARNING:SIM106", "Parameter {0} corresponds both as a non zero parameter \
                            and a rate rule, verify behavior".format(element))
                            removeParameters.append(element)
            # it is an assigment rule
            elif rawArule[2] is True:
                '''
                 normal species observables references in functions keep the format <speciesName>_<compartment> in function references,
                 and observables dict keeps track of that. however when a species is defined by an assignment function we wish to 
                 keep track of reference <speciesName> that points to a standard BNGL function
                '''
                #if rawArule[0] in observablesDict:
                #    del observablesDict[rawArule[0]]

                # it was originially defined as a zero parameter, so delete it from the parameter list definition                
                if rawArule[0] in zRules:
                    # dont show assignment rules as parameters
                    zRules.remove(rawArule[0])
                    #zRules.append([rawArule[0] + '_assignment', rawArule[1], rawArule[2], rawArule[3]])

                    #aParameters[rawArule[0]] = 'arj' + rawArule[0]
                    #tmp = list(rawArule)
                    #tmp[0] = 'arj' + rawArule[0]
                    #rawArule= tmp
                    matches = [molecules[x] for x in molecules if molecules[x]['name'] == rawArule[0]]
                    if matches:
                        if matches[0]['isBoundary']:
                            artificialObservables[rawArule[0] + '_ar'] = writer.bnglFunction(rawArule[1][0],rawArule[0]+'_ar()',[],compartments=compartmentList,reactionDict=self.reactionDictionary)
                            continue
                        else:
                            logMess('ERROR:SIM201', 'Variables that are both changed by an assignment rule and reactions are not \
                            supported in BioNetGen simulator. The variable will be split into two'.format(rawArule[0]))
                            artificialObservables[rawArule[0] + '_ar'] = writer.bnglFunction(rawArule[1][0],rawArule[0]+'_ar()',[],compartments=compartmentList,reactionDict=self.reactionDictionary)
                            continue
                    elif rawArule[0] in [observablesDict[x] for x in observablesDict]:
                        artificialObservables[rawArule[0] + '_ar'] = writer.bnglFunction(rawArule[1][0],rawArule[0]+'_ar()',[],compartments=compartmentList,reactionDict=self.reactionDictionary)
                        continue

                elif rawArule[0] in molecules:
                    if molecules[rawArule[0]]['isBoundary']:
                        artificialObservables[rawArule[0]+'_ar'] = writer.bnglFunction(rawArule[1][0],rawArule[0]+'_ar()',[],compartments=compartmentList,reactionDict=self.reactionDictionary)
                        continue
                else:
                    #check if it is defined as an observable
                    candidates =  [idx for idx,x in enumerate(observablesDict) if rawArule[0] == x]
                    assigObsFlag = False
                    for idx in candidates:
                        #if re.search('\s{0}\s'.format(rawArule[0]),observables[idx]):
                        artificialObservables[rawArule[0]+ '_ar'] = writer.bnglFunction(rawArule[1][0],rawArule[0]+'_ar()',[],compartments=compartmentList,reactionDict=self.reactionDictionary)
                        assigObsFlag = True
                        break
                    if assigObsFlag:
                        continue
                # if its not a param/species/observable
                artificialObservables[rawArule[0]] = writer.bnglFunction(rawArule[1][0],rawArule[0]+'()',[],compartments=compartmentList,reactionDict=self.reactionDictionary)

            else:
                '''
                if for whatever reason you have a rule that is not assigment
                or rate and it is initialized as a non zero parameter, give it
                a new name
                '''
                if rawArule[0] not in zparams:
                    ruleName = 'ar' + rawArule[0]
                else:
                    ruleName = rawArule[0]
                    zRules.remove(rawArule[0])
                arules.append(writer.bnglFunction(rawArule[1][0], ruleName, [], compartments=compartmentList, reactionDict=self.reactionDictionary))
                aParameters[rawArule[0]] = 'ar' + rawArule[0]
            '''
            elif rawArule[2] == True:
                for parameter in parameters:
                    if re.search('^{0}\s'.format(rawArule[0]),parameter):
                        print '////',rawArule[0]
            '''
            #arules.append('%s = %s' %(rawArule[0],newRule))
        return aParameters, arules, zRules, artificialReactions, removeParameters, artificialObservables

    def convertToStandardUnits(self, parameterValue, unitDefinition):

        for factor in unitDefinition:
            if factor['scale'] != 0:
                parameterValue *= 10 ** factor['scale']
            if factor['multiplier'] !=  1:
                parameterValue *= factor['multiplier']
            if factor['exponent'] != 1:
                parameterValue **= factor['exponent']
        return parameterValue


    def convertToStandardUnitString(self, parameterValue, unitDefinition):
        for factor in unitDefinition:
            if factor['multiplier'] !=  1:
                parameterValue = '({0} * {1})'.format(parameterValue, factor['multiplier'])
            if factor['exponent'] != 1:
                parameterValue = '({0} ^ {1})'.format(parameterValue, factor['exponent'])
            if factor['scale'] != 0:
                parameterValue = '({0} * 1e{1})'.format(parameterValue, factor['scale'])
        return parameterValue



    def __getRawParameters(self, parameter):
        parameterSpecs = {}
        parameterSpecs['id'] = parameter.getId()
        parameterSpecs['value'] = parameter.getValue()
        parameterSpecs['name'] = parameter.getName()
        parameterSpecs['units'] = parameter.getUnits()

        return parameterSpecs

    def getParameters(self):
        parameters = []
        zparam = []
        for parameter in self.model.getListOfParameters():
            parameterSpecs = (parameter.getId(), parameter.getValue(), parameter.getConstant(), parameter.getUnits())
            #print parameterSpecs
            #reserved keywords
            if parameterSpecs[0] == 'e':
                parameterSpecs = ('are', parameterSpecs[1])
            if parameterSpecs[1] == 0:
                zparam.append(parameterSpecs[0])
            else:
                """
                if self.convertSubstanceUnits:
                    newParameterSpecs = [parameterSpecs[0], self.convertToStandardUnits(parameterSpecs[1], self.unitDictionary[parameterSpecs[3]]), parameterSpecs[2], parameterSpecs[3]]
                    parameters.append('{0} {1} #original units:{2}={3}'.format(newParameterSpecs[0], newParameterSpecs[1], newParameterSpecs[3],parameterSpecs[1]))

                else:
                """
                if parameter.getUnits() != '':
                    parameters.append('{0} {1} #units:{2}'.format(parameterSpecs[0], parameterSpecs[1], parameter.getUnits()))
                else:
                    parameters.append('{0} {1}'.format(parameterSpecs[0], parameterSpecs[1]))

        #return ['%s %f' %(parameter.getId(),parameter.getValue()) for parameter in self.model.getListOfParameters() if parameter.getValue() != 0], [x.getId() for x in self.model.getListOfParameters() if x.getValue() == 0]
        return parameters, zparam

    def getSpecies(self, translator={}, parameters=[]):
        '''
        in sbml parameters and species have their own namespace. not so in
        bionetgen, so we need to rename things if they share the same name
        '''
        def default_to_regular(d):
            if isinstance(d, defaultdict):
                d = {k: default_to_regular(v) for k, v in d.iteritems()}
            return d

        #find concentration units
        unitDefinitions = self.getUnitDefinitions()

        if 'substance' in unitDefinitions:
            substance = '*'.join([x['name'] for x in unitDefinitions['substance']])
        else:
            substance = 'mol'
        if 'volume' in unitDefinitions:
            volume = '/'.join([x['name'] for x in unitDefinitions['volume']])
        else:
            volume = 'L'
        concentrationUnits = '{0}/{1}'.format(substance,volume)
        annotationMoleculesText = {}
        moleculesText = []
        speciesText = []
        observablesText = []
        observablesDict = {}
        names = []
        rawSpeciesName = translator.keys()
        speciesTranslationDict = {}
        compartmentDict = {}
        compartmentDict[''] = 1
        speciesAnnotationInfo = default_to_regular(self.getFullAnnotation())
        annotationInfo = {'moleculeTypes': {}, 'species': {}}
        for compartment in self.model.getListOfCompartments():
            compartmentDict[compartment.getId()] = compartment.getSize()
        #speciesText.append('#units: {0}'.format(concentrationUnits))
        unitFlag = True
        for species in self.model.getListOfSpecies():
            rawSpecies = self.getRawSpecies(species, parameters)
            #if rawSpecies['returnID'] in self.boundaryConditionVariables:
            #    continue
            if (rawSpecies['compartment'] != ''):
                self.tags[rawSpecies['identifier']] = '@%s' % (rawSpecies['compartment'])
            if(rawSpecies['returnID'] in translator):
                if rawSpecies['returnID'] in rawSpeciesName:
                    rawSpeciesName.remove(rawSpecies['returnID'])
                if translator[rawSpecies['returnID']].getSize() == 1 \
                    and translator[rawSpecies['returnID']].molecules[0].name not in names \
                        and translator[rawSpecies['returnID']].molecules[0].name not in rawSpeciesName:
                    names.append(translator[rawSpecies['returnID']].molecules[0].name)
                    annotationTemp = []
                    if rawSpecies['returnID'] in speciesAnnotationInfo:
                        for annotation in speciesAnnotationInfo[rawSpecies['returnID']]:
                            parts = annotation.split('_')
                            header = annotationHeader[parts[0]]
                            qual = parts[1].lower() + ''.join([x.capitalize() for x in parts[2:]])
                            entry = ', '.join([':'.join(x.split('/')[-2:]) for x in speciesAnnotationInfo[rawSpecies['returnID']][annotation]])
                            annotationTemp.append('#^ {0}:{1} {2}'.format(header,qual, entry))
                    
                    moleculesText.append(translator[rawSpecies['returnID']].str2())
                    if rawSpecies['returnID'] in speciesAnnotationInfo:
                        annotationInfo['moleculeTypes'][translator[rawSpecies['returnID']].str2()] = annotationTemp
                        del speciesAnnotationInfo[rawSpecies['returnID']]
            else:
                moleculesText.append(rawSpecies['returnID'] + '()')
                if rawSpecies['returnID'] in speciesAnnotationInfo:
                    annotationInfo['moleculeTypes'][rawSpecies['returnID']] = speciesAnnotationInfo[rawSpecies['returnID']]
                    del speciesAnnotationInfo[rawSpecies['returnID']]

            temp = '$' if rawSpecies['isConstant'] != 0 else ''
            tmp = translator[str(rawSpecies['returnID'])] if rawSpecies['returnID'] in translator \
                else rawSpecies['returnID'] + '()'
            if rawSpecies['initialConcentration'] >= 0 or rawSpecies['initialAmount'] >=0:
                tmp2 = temp
                if rawSpecies['identifier'] in self.tags:
                    tmp2 = (self.tags[rawSpecies['identifier']])
                if rawSpecies['initialAmount'] > 0.0:
                    speciesText.append('{0}:{1}{2} {3} #{4} #{5}'.format(tmp2, temp, str(tmp), rawSpecies['initialAmount'],rawSpecies['returnID'],rawSpecies['identifier']))
                elif rawSpecies['initialConcentration'] > 0.0:
                    if self.isConversion:
                        # convert to molecule counts
                        if 'substance' in unitDefinitions:
                            newParameterStr = self.convertToStandardUnitString(rawSpecies['initialConcentration'], unitDefinitions['substance'])
                            newParameter = self.convertToStandardUnits(rawSpecies['initialConcentration'], unitDefinitions['substance']) #conversion to moles
                        else:
                            newParameter = rawSpecies['initialConcentration']
                            newParameterStr = str(rawSpecies['initialConcentration'])
                        newParameter = newParameter * 6.022e23 # convertion to molecule counts                        
                        #get compartment size
                        compartmentSize = self.model.getCompartment(rawSpecies['compartment']).getSize()
                        newParameter = compartmentSize * newParameter
                        if unitFlag:
                            speciesText.append('{0}:{1}{2} {3} # {4}mol/L * 6.022e23/mol *{7}L #{5} #{6}'.format(tmp2, temp, str(tmp), 
                                                                                    newParameter,newParameterStr,rawSpecies['returnID'],
                                                                                    rawSpecies['identifier'], compartmentSize, concentrationUnits))
                            unitFlag = False
                        else:
                            speciesText.append('{0}:{1}{2} {3} #original {4}{8}  #{5} #{6}'.format(tmp2, temp, str(tmp), 
                                                                                    newParameter,rawSpecies['initialConcentration'],rawSpecies['returnID'],
                                                                                    rawSpecies['identifier'], compartmentSize, concentrationUnits))
                    else:
                        speciesText.append('{0}:{1}{2} {3} #{4} #{5}'.format(tmp2, temp, str(tmp), rawSpecies['initialConcentration'],rawSpecies['returnID'],rawSpecies['identifier']))

                elif rawSpecies['isConstant']:
                    speciesText.append('{0}:{1}{2} {3} #{4} #{5}'.format(tmp2, temp, str(tmp), 0,rawSpecies['returnID'],rawSpecies['identifier']))
            if rawSpecies['returnID'] == 'e':
                modifiedName = 'are'
            else:
                modifiedName = rawSpecies['returnID']
            # user defined zero molecuels are not included in the observable list
            if str(tmp) != '0':
                if rawSpecies['compartment'] != '' and len(list(self.model.getListOfCompartments())) > 1:
                    observablesText.append('Species {0}_{3} @{3}:{1} #{2}'.format(modifiedName, tmp,rawSpecies['name'],rawSpecies['compartment']))
                else:
                    observablesText.append('Species {0}_{3} @{3}:{1} #{2}'.format(modifiedName, tmp,rawSpecies['name'],rawSpecies['compartment']))
                observablesDict[modifiedName] = '{0}_{1}'.format(modifiedName,rawSpecies['compartment'])
                speciesTranslationDict[rawSpecies['identifier']] = tmp
        sorted(rawSpeciesName,key=len)
        for species in rawSpeciesName:
            if translator[species].getSize()==1 and translator[species].molecules[0].name not in names:
                names.append(translator[species].molecules[0].name)
                moleculesText.append(translator[species].str2())

        annotationInfo['species'] = speciesAnnotationInfo

        #moleculesText.append('NullSpecies()')
        #speciesText.append('$NullSpecies() 1')

        self.speciesMemory = []
        return list(set(moleculesText)),speciesText,observablesText,speciesTranslationDict, observablesDict, annotationInfo

    def getInitialAssignments(self, translator, param, zparam, molecules, initialConditions):
        '''
        process the get initial assignments section. This can be used to initialize
        parameters or species, so we have to account for both checking both arrays
        '''
        param2 = param
        zparam2 = zparam
        initialConditions2 = initialConditions
        pparam = {}
        for element in param:
            pparam[element.split(' ')[0]] =(element.split(' ')[1],None)
        for element in zparam:
            pparam[element] = (0,None)
        for species in self.model.getListOfSpecies():
            tmp = self.getRawSpecies(species)

            #name = species.getName() if species.isSetName() else species.getId()
            name = tmp['returnID']
            constant = '$' if species.getConstant() or species.getBoundaryCondition() else ''
            if name in  translator:
                extendedStr = '@{0}:{2}{1}'.format(species.getCompartment(),translator[name],constant)
            else:
                extendedStr = '@{0}:{2}{1}()'.format(tmp['compartment'],standardizeName(tmp['name']),constant)
            initConc = species.getInitialConcentration() if \
            species.isSetInitialConcentration() else species.getInitialAmount()
            pparam[species.getId()] = (initConc,extendedStr)
        from copy import copy
        for initialAssignment in self.model.getListOfInitialAssignments():
            symbol = initialAssignment.getSymbol()
            math = libsbml.formulaToString(initialAssignment.getMath())
            for element in pparam:
                if element in math:
                    math = re.sub(r'(\W|^)({0})(\W|$)'.format(element),
                    r'\1({0})\3'.format(pparam[element][0]),math)

            # removing non bngl math elements  for their equivalents
            math =  writer.bnglFunction(math,'',[]).split(' = ')[1]
            param2 = [x for x in param if '{0} '.format(symbol) not in x]

            zparam2 = [x for x in zparam if symbol != x]
            '''
            if (len(param2) != len(param)) or (len(zparam) != len(zparam2)):
                param2.append('{0} {1}'.format(symbol,math))
                param = param2
                zparam = zparam2
            '''
            try:
                if pparam[symbol][1] == None:
                    param2.append('{0} {1}'.format(symbol,math))
                    param = param2
                    zparam = zparam2
                else:
                    
                    initialConditions2 = [x for x in initialConditions if '#{0}'.format(symbol) not in x]
                    initialConditions2.append('{0} {1} #{2}'.format(pparam[symbol][1],math,symbol))
                    initialConditions = initialConditions2
            except:
                continue

        return param,zparam,initialConditions
            
            
    def getSpeciesAnnotation(self):
        if self.speciesAnnotation:
            return self.speciesAnnotation
        
        self.speciesAnnotation = defaultdict(list)

        for species in self.model.getListOfSpecies():
            rawSpecies = self.getRawSpecies(species,logEntries=False)
            annotationXML = species.getAnnotation()
            lista = libsbml.CVTermList()
            libsbml.RDFAnnotationParser.parseRDFAnnotation(annotationXML,lista)
            if lista.getSize() == 0:
                self.speciesAnnotation[rawSpecies['returnID']] =  []
            else:
                for idx in range(lista.getSize()):
                    self.speciesAnnotation[rawSpecies['returnID']].append(lista.get(idx).getResources())
       
        return self.speciesAnnotation

    def getFullAnnotation(self):
        if self.speciesAnnotationDict:
            return self.speciesAnnotationDict
        self.speciesAnnotationDict = defaultdict(lambda: defaultdict(list))

        for species in self.model.getListOfSpecies():
            rawSpecies = self.getRawSpecies(species, logEntries=False)
            annotationXML = species.getAnnotation()
            lista = libsbml.CVTermList()
            libsbml.RDFAnnotationParser.parseRDFAnnotation(annotationXML, lista)
            if lista.getSize() == 0:
                continue
            else:
                for idx in range(lista.getSize()):
                    for idx2 in range(0, lista.get(idx).getResources().getLength()):
                        resource = lista.get(idx).getResources().getValue(idx2)
                        qualifierType = lista.get(idx).getQualifierType()
                        qualifierDescription = bioqual[lista.get(idx).getBiologicalQualifierType()] if qualifierType \
                            else modqual[lista.get(idx).getModelQualifierType()]
                        self.speciesAnnotationDict[rawSpecies['returnID']][qualifierDescription].append(resource)

        return self.speciesAnnotationDict

    def getModelAnnotation(self):
        modelAnnotation = []
        annotationXML = self.model.getAnnotation()
        lista = libsbml.CVTermList()
        libsbml.RDFAnnotationParser.parseRDFAnnotation(annotationXML, lista)
        if lista.getSize() == 0:
            modelAnnotations = []
        else:
            tempDict = {}
            for element in [2, 3, 4, 5, 6]:
                if lista.get(element) == None:
                    continue
                tempDict[lista.get(element).getBiologicalQualifierType()] = lista.get(element)
            if 3 in tempDict:
                modelAnnotation = tempDict[3].getResources()
            elif 0 in tempDict and ('GO' in tempDict[0].getResources().getValue(1) or 'kegg' in tempDict[0].getResources().getValue(1)):
                modelAnnotation = tempDict[0].getResources()
            elif 5 in tempDict:
                modelAnnotation = tempDict[5].getResources()
            else:
                if lista.get(3) != None and ('GO' in lista.get(3).getResources().getValue(0) or 'kegg' in lista.get(3).getResources().getValue(0)):
                    modelAnnotation = lista.get(3).getResources()
                    
                elif lista.get(4) != None and ('GO' in lista.get(4).getResources().getValue(0) or 'kegg' in lista.get(4).getResources().getValue(0)):
                    modelAnnotation = lista.get(4).getResources()
                elif lista.get(5) != None and ('GO' in lista.get(5).getResources().getValue(0) or 'kegg' in lista.get(5).getResources().getValue(0)):
                    modelAnnotation = lista.get(5).getResources()
                else:
                    if lista.get(3) != None and ('reactome' in lista.get(3).getResources().getValue(0)):
                        modelAnnotation = lista.get(3).getResources()
                        
                    elif lista.get(4) != None and ('reactome' in lista.get(4).getResources().getValue(0)):
                        modelAnnotation = lista.get(4).getResources()
                    elif lista.get(5) != None and ('reactome' in lista.get(5).getResources().getValue(0)):
                        modelAnnotation = lista.get(5).getResources()
                    elif lista.get(2) != None:
                        modelAnnotation = lista.get(2).getResources()
        return modelAnnotation
        
    def getSpeciesInfo(self, name):
        return self.getRawSpecies(self.model.getSpecies(name))

    def writeLog(self, translator):
        rawSpecies = [self.getRawSpecies(x) for x in self.model.getListOfSpecies()]
        log['species'].extend([x[0] for x in rawSpecies if x[0] not in translator])
        logString = ''
        #species stuff
        if(len(log['species']) > 0):
            logString += "Species we couldn't recognize:\n"
            for element in log['species']:
                logString += '\t%s\n' % element
        if(len(log['reactions']) > 0):
            logString += "Reactions we couldn't infer more about due to \
            insufficient information:"
            for element in log['reactions']:
                logString += '\t%s + %s -> %s\n' % (element[0][0],
                                                    element[0][1],
                                                    element[1])
        return logString

    def getStandardName(self,name):
        if name in self.speciesDictionary:
            return self.speciesDictionary[name]
        return name
    
    
def standardizeName(name):
    '''
    Remove stuff not used by bngl
    '''
    name2 = name
    
    sbml2BnglTranslationDict = {"^":"",
                                "'":"",
                                "*":"m"," ":"_",
                                "#":"sh",
                                ":":"_",'α':'a',
                                'β':'b',
                                'γ':'g',"(":"__",
                                ")":"__",
                                " ":"","+":"pl",
                                "/":"_",":":"_",
                                "-":"_",
                                ".":"_",
                                '?':"unkn",
                                ',':'_',
                                '[':'__',
                                  ']':'__',
                                  '>':'_',
                                  '<':'_'}
                                
    for element in sbml2BnglTranslationDict:
        name = name.replace(element,sbml2BnglTranslationDict[element])
    name = re.sub('[\W]', '', name)
    return name
