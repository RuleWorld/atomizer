import networkx as nx
import stateTransitionDiagram as std
import argparse


def defineConsole():
    """
    defines the program console line commands
    """
    parser = argparse.ArgumentParser(description='SBML to BNGL translator')
    parser.add_argument('-i', '--input', type=str, help='settings file', required=True)
    return parser


def getCounter():
    if not hasattr(getCounter, 'counter'):
        getCounter.counter = 0

    getCounter.counter += 1
    return getCounter.counter


def createNode(graph, name, graphicsDict, labelGraphicsDict, isGroup, gid):
    idNumber = getCounter()
    if isGroup:
        graph.add_node(name, graphics=graphicsDict, LabelGraphics=labelGraphicsDict, isGroup=isGroup, id=idNumber)
    else:
        graph.add_node(name, graphics=graphicsDict, LabelGraphics=labelGraphicsDict, gid=gid, id=idNumber)


def generateSTD(nodeList, edgeList):
    graph = nx.DiGraph()

    for molecule in nodeList:
        createNode(graph, molecule, {'type': 'roundrectangle', 'fill': '#FFDD99'},
                   {'fontSize': 16, 'fontStyle': "bold", 'alignment': "right", 'autoSizePolicy': "node_size"}, 1, 0)

        for node in nodeList[molecule]:
            componentLegend = '/'.join([x[0] for x in node])

            nodeName = u''
            nodeId = []
            for bit in node:
                if bit[1]:
                    nodeName += u"\u25CF"
                    nodeId.append(bit[0])
                else:
                    nodeName += u"\u25CB"
            createNode(graph, '{0}_{1}'.format(molecule, '/'.join(nodeId)), {'type': "roundrectangle"}, {'text': nodeName}, 0, graph.node[molecule]['id'])

        createNode(graph, '{0}_legend'.format(molecule), {}, {'text': componentLegend}, 0, graph.node[molecule]['id'])
 
    for molecule in nodeList:
        for edge in edgeList[molecule]:
            # print edge
            nodeId0 = [x[0] for x in edge[0] if x[1]]
            nodeId1 = [x[0] for x in edge[1] if x[1]]
            bidirectional = (edge[1], edge[0]) in edgeList[molecule]
            if ('{0}_{1}'.format(molecule, '/'.join(nodeId1)), '{0}_{1}'.format(molecule, '/'.join(nodeId0))) not in graph.edges():
                if bidirectional:
                    graph.add_edge('{0}_{1}'.format(molecule, '/'.join(nodeId0)), '{0}_{1}'.format(molecule, '/'.join(nodeId1)),
                                   graphics={'fill': '#000000', 'sourceArrow': "standard", 'targetArrow': "standard"})

                else:
                    graph.add_edge('{0}_{1}'.format(molecule, '/'.join(nodeId0)), '{0}_{1}'.format(molecule, '/'.join(nodeId1)),
                                   graphics={'fill': '#000000', 'targetArrow': "standard"})

    return graph


def outputGraph(graph, fileName):
    nx.write_gml(graph, fileName)

if __name__ == "__main__":
    parser = defineConsole()
    namespace = parser.parse_args()
    inputFile = namespace.input

    nodeList, edgeList = std.getContextRequirements(inputFile)
    graph = generateSTD(nodeList, edgeList)
    outputGraph(graph, 'output.gml')
