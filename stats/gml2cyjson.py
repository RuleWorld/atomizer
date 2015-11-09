import random
from copy import copy

def processEdge(edgeInfo):
    pass

def gml2cyjson(gmlText, graphtype=None):
    """
    Converts a gml graph definition to the format that cytoscape.js expects
    """
    r = lambda: random.randint(0, 255)

    jsonDict = {}
    jsonDict['style'] =  [
    {
      'selector': 'node',
      'css': {
        'content': 'data(label)',
        'text-valign': 'center',
        'text-halign': 'center'
      }
    },
    {
      'selector': '$node > node',
      'css': {
        'padding-top': '20px',
        'padding-left': '10px',
        'padding-bottom': '10px',
        'padding-right': '10px',
        'text-valign': 'top',
        'text-halign': 'center'
      }
    },
    {
      'selector': 'edge',
      'css': {
        'target-arrow-shape': 'none'
      }
    },
    {
      'selector': ':selected',
      'css': {
        'background-color': 'black',
        'line-color': 'black',
        'target-arrow-color': 'black',
        'source-arrow-color': 'black'
      }
    }
  ]
    jsonDict['elements'] = {'nodes':[]}
    jsonDict['elements']['edges'] = []
    #for nd in gmlText.node:
    #    gmlText.node[nd].pop('id')
    #jsgrph = json_graph.node_link_data(gmlText)
    colorDict = {}
    for node in gmlText.node:
        if gmlText.node[node] == {}:
            continue
        tmp = {'data':{}}
        tmp['data']['id'] = str(node)
        tmp['data']['label'] = str(gmlText.node[node]['label'])
        
        if 'gid' in gmlText.node[node]:
            tmp['data']['parent'] =  str(gmlText.node[node]['gid'])
            if str(gmlText.node[node]['gid']) not in colorDict:
                if 'gid' in gmlText.node[str(gmlText.node[node]['gid'])]:
                    if str(gmlText.node[str(gmlText.node[node]['gid'])]['gid']) not in colorDict:
                        if graphtype == 'regulatory':
                            newColor = gmlText.node[node]['graphics']['fill']
                        else:
                            newColor = '#%02X%02X%02X' % (r(), r(), r())
                        colorDict[str(gmlText.node[str(gmlText.node[node]['gid'])]['gid'])] = newColor
                        colorDict[str(gmlText.node[node]['gid'])] = newColor
                    else:
                        colorDict[str(gmlText.node[node]['gid'])] = colorDict[str(gmlText.node[str(gmlText.node[node]['gid'])]['gid'])]
                else:
                    if graphtype == 'regulatory':
                        colorDict[str(gmlText.node[node]['gid'])] = gmlText.node[node]['graphics']['fill']
                    else:
                        colorDict[str(gmlText.node[node]['gid'])] = '#%02X%02X%02X' % (r(), r(), r())
            colorDict[str(node)] = colorDict[str(gmlText.node[node]['gid'])]
        if str(node) not in colorDict:
            if graphtype == 'regulatory':
                colorDict[str(node)] = gmlText.node[node]['graphics']['fill']
            else:
                colorDict[str(node)] = '#%02X%02X%02X' % (r(), r(), r())
        tmp['data']['faveColor'] = colorDict[str(node)]

        jsonDict['elements']['nodes'].append(tmp)
    for link in gmlText.edge:
        for dlink in gmlText.edge[link]:
            if link != '' and dlink != '':
                tmp = {'data':{}}
                tmp['data']['source'] = int(link)
                tmp['data']['target'] = int(dlink)
                tmp['data']['id'] = '{0}_{1}'.format(tmp['data']['source'], tmp['data']['target'])


                if graphtype == 'regulatory':
                    if 'graphics' in gmlText.edge[link][dlink]:
                        tmp['data']['faveColor'] = gmlText.edge[link][dlink]['graphics']['fill']
                        jsonDict['elements']['edges'].append(tmp)
                    else:
                        for multiedge in gmlText.edge[link][dlink]:
                            if 'graphics' in gmlText.edge[link][dlink][multiedge]:
                                tmp['data']['faveColor'] = gmlText.edge[link][dlink][multiedge]['graphics']['fill']
                                jsonDict['elements']['edges'].append(copy(tmp))
                else:
                    tmp['data']['faveColor'] = colorDict[str(link)]
                    jsonDict['elements']['edges'].append(tmp)

    jsonDict['layout'] = {
    'name': 'cose',
    'padding': 4,
     'fit'                 : True, 
   'nodeRepulsion'       : 10000,
    'nodeOverlap'         : 10,
    'idealEdgeLength'     : 10,
   'edgeElasticity'      : 100,
    'nestingFactor'       : 5, 
    'gravity'             : 250, 
    'numIter'             : 100,
    'initialTemp '        : 200,
    'coolingFactor'       : 0.95, 
    'minTemp'             : 1  },
  
    jsonDict['ready'] =  'function(){window.cy = this;}'

    return jsonDict

'''
if __name__ == '__main__':
    import networkx as nx
    import pprint
    #with open ('/home/proto/workspace/bionetgen/bng2/Models2/toy-jim_regulatory.gml','r') as f:
    #    s = f.read()
    s = nx.read_gml('/home/proto/workspace/bionetgen/bng2/Models2/toy-jim_regulatory.gml')    
    graph = gml2cyjson(s,'regulatory')
    pprint.pprint(graph['elements']['edges'])
''' 

