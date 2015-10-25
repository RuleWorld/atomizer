from utils import readBNGXML
import argparse
import collections
import csv
import stateTransitionDiagram as std
import pprint
import xlwt
import yaml

class PrettyDefaultDict(collections.defaultdict):
    __repr__ = dict.__repr__

def getParameterDictionary(bnglNamespace):

    parameterDict = PrettyDefaultDict(lambda: PrettyDefaultDict(set))
    simpleParameterDict = collections.defaultdict(set)

    labels, centers, contexts, products, atomicArrays, actions, doubleAction = std.extractCenterContext(bnglNamespace['rules'])

    for label, center, product, rule in zip(labels, centers, products, bnglNamespace['rules']):
        for pattern in product:
            for element in pattern:
                for molecule in element.split('.'):
                    for rate in rule[0].rates:
                        value = bnglNamespace['parameters'][rate] if rate in bnglNamespace['parameters'] else 'nma'
                        moleculeName = molecule.split('(')[0].split('%')[0]
                        component = molecule.split('(')[1].strip(')')
                        parameterDict[moleculeName][component].add(value)
                        simpleParameterDict[molecule].add(value)

    return parameterDict

def defineConsole():
    """
    defines the program console line commands
    """
    parser = argparse.ArgumentParser(description='SBML to BNGL translator')
    parser.add_argument('-i', '--input', type=str, help='settings file', required=True)
    return parser


if __name__ == "__main__":
    parser = defineConsole()
    namespace = parser.parse_args()
    inputFile = namespace.input
    modelNameList = ['egfr/output19.xml', 'egfr/output48.xml', 'egfr/output151.xml', 'egfr/output543.xml']
    #modelName = ['egfr/output151.xml']
    parameterSpace = []
    for element in modelNameList:
        parameterSpace.append(getParameterDictionary(readBNGXML.parseFullXML(element)))

    keys = collections.Counter([y for x in parameterSpace for y in x])
    wb = xlwt.Workbook()

    ws = wb.add_sheet('Units')
    unitColumn = []
    for midx, element in enumerate(modelNameList):
        modelName = element.split('.')[0]
        ymlName = modelName + '.bngl.yml'

        try:
            with open(ymlName, 'r') as f:
                annotationDict = yaml.load(f)
        except IOError:
            continue
        ws.write(midx + 1, 0, modelName)
        for unit in annotationDict['units']:
            if unit not in unitColumn:
                unitColumn.append(unit)
                ws.write(0, unitColumn.index(unit) + 1, unit)
            #for instance in annotationDict['units'][unit]:
            instance = annotationDict['units'][unit][0]
            if instance['exponent'] == 1:
                ws.write(midx + 1, unitColumn.index(unit) + 1, '{0} ({1})'.format(instance['multiplier'], instance['name']))
            else:
                ws.write(midx + 1, unitColumn.index(unit) + 1, '{0}^{2} ({1})'.format(instance['multiplier'], instance['name'], instance['exponent']))
    for molecule in keys:
        if keys[molecule] == 1:
            continue
        ws = wb.add_sheet(molecule)
        componentColumn = []
        for midx, model in enumerate(parameterSpace):
            ws.write(midx + 1, 0, modelNameList[midx])
            for component in model[molecule]:
                if component not in componentColumn:
                    componentColumn.append(component)
                    if '!' in component:
                        componentTxt = component.split('!')[0] + ('(association)')
                    elif '~0' in component:
                        componentTxt = component.split('~')[0] + ('(repression)')
                    elif '~' in component:
                        componentTxt = component.split('~')[0] + ('(activation)')
                    else:
                        componentTxt = component + ('(dissociation)')

                    ws.write(0, componentColumn.index(component) + 1, componentTxt)
                data = '/'.join(model[molecule][component]) 
                ws.write(midx+1, componentColumn.index(component) + 1, data)

    wb.save('19_151_48.xls')
    '''
    with open('19_151_48.csv','wb') as f:
        writer = csv.DictWriter(f, keys)
        writer.writeheader()
        writer.writerows(parameterSpace)
    '''
    