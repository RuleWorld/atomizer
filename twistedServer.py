                # -*- coding: utf-8 -*-
"""
Created on Fri May 2 16:56:13 2014

@author: proto
"""

import os

import sys
# Restrict to a particular path.
from twisted.web import xmlrpc, server
from twisted.internet import reactor, task
import threading
import SBMLparser.utils.consoleCommands as consoleCommands

import tempfile
import gml2sbgn.libsbgn as libsbgn
import networkx
from subprocess import call
import yaml
from stats.gml2cyjson import gml2cyjson
import json
import StringIO

sys.path.insert(0, 'SBMLparser')
import SBMLparser.libsbml2bngl as libsbml2bngl
import SBMLparser.utils.readBNGXML as readBNGXML
import SBMLparser.utils.annotationExtender as annotationExtender
import SBMLparser.utils.nameNormalizer as normalizer
import SBMLparser.utils.modelComparison as modelComparison
import SBMLparser.rulifier.stdgraph as stdgraph
bngDistro = '/home/ubuntu/wokspace/bionetgen/bng2/BNG2.pl'
iid = 1
iid_lock = threading.Lock()

#bngDistro = '/home/proto/workspace/bionetgen/bng2/BNG2.pl'


def next_id():
    global iid
    with iid_lock:
        result = iid
        iid += 1
    return result

processDict = {}



def freeQueue(ticket):
    processDict.pop(ticket)

class AtomizerServer(xmlrpc.XMLRPC):

    def addToDict(self, ticket, result):
        if len(processDict) > 40:
            s = min(processDict.keys())
            processDict.pop(s)
        processDict[ticket] = result

    def atomize(self, ticket, xmlFile, atomize, userConf = None):
        reaction = 'config/reactionDefinitions.json'
        try:
            logStream = StringIO.StringIO()
            if userConf:
                jsonpointer = tempfile.mkstemp(suffix='.json', text=True)
                with open(jsonpointer[1], 'w') as f:
                    f.write(userConf)
                jsonpointer = jsonpointer[1]
            else:
                jsonpointer = None
            result = libsbml2bngl.readFromString(xmlFile,
                                                 reaction, False, jsonpointer, atomize, logStream)

            if result and atomize:
                pointer = tempfile.mkstemp(suffix='.bngl', text=True)
                with open(pointer[1], 'w') as f:
                    f.write(result.finalString)
                print pointer[1]
                bnglresult = libsbml2bngl.postAnalyzeString(pointer[1], bngDistro, result.database)
            else:
                bnglresult = result.finalString
            self.addToDict(ticket, [bnglresult, logStream.getvalue(), {'finalspecies':result.database.species, 'rawreactions':result.database.rawreactions}])
            print 'success', ticket

        except:
            self.addToDict(ticket, -5)
            print 'failure', ticket
        finally:
            task.deferLater(reactor, 600,  freeQueue, ticket)


    def extractMoleculeTypes(self,ticket,bnglContents, bnglContents2):

        moleculeTypesList = []
        for element in [bnglContents, bnglContents2]:
            pointer = tempfile.mkstemp(suffix='.bngl', text=True)
            with open(pointer[1], 'w') as f:
                f.write(element)
            print pointer[1]
            consoleCommands.setBngExecutable(bngDistro)
            consoleCommands.bngl2xml(pointer[1])

            xmlFileName = pointer[1].split('.')[0] + '.xml'
            xmlFileName = xmlFileName.split(os.sep)[-1]
            moleculeTypes, _, _ = readBNGXML.parseXML(xmlFileName)
            moleculeTypesList.append(moleculeTypes)
            os.remove(xmlFileName)
        self.addToDict(ticket, moleculeTypesList)
        print 'success', ticket

    def compareFiles(self, ticket, bnglContents, bnglContents2, mappingFile):
        finalBNGLContent = []
        finalNamespace = []
        try:
            for mapInfo, bnglContent in zip(mappingFile['model'], [bnglContents, bnglContents2]):
                pointer = tempfile.mkstemp(suffix='.bngl', text=True)
                with open(pointer[1], 'w') as f:
                    f.write(bnglContent)

                print pointer[1]
                consoleCommands.setBngExecutable(bngDistro)
                consoleCommands.bngl2xml(pointer[1])

                xmlFileName = pointer[1].split('.')[0] + '.xml'
                xmlFileName = xmlFileName.split(os.sep)[-1]
                bnglNamespace = readBNGXML.parseFullXML(xmlFileName)
                normalizer.normalizeNamespace(bnglNamespace, mapInfo)
                
                finalBNGLContent.append(readBNGXML.createBNGLFromDescription(bnglNamespace))
                finalNamespace.append(bnglNamespace)



                # os.remove(pointer[1])
                os.remove(xmlFileName)

            similarity = modelComparison.evaluateSimilarity(finalNamespace[0], finalNamespace[1])
            self.addToDict(ticket, [finalBNGLContent, similarity])
            print 'success', ticket
        except:
            self.addToDict(ticket,-5)
            print 'failure',ticket
        finally:
            task.deferLater(reactor, 600,  freeQueue, ticket)
            


        pass


    def generateAnnotation(self, ticket, xmlFile):

        print ticket
        reaction = 'config/reactionDefinitions.json'

        pointer = tempfile.mkstemp(suffix='.xml', text=True)
        with open(pointer[1], 'w') as f:
            f.write(xmlFile)
        '''
            call(['python','annotationExtender.py',
            '-i',pointer[1],
            '-o',pointer[1]+'.xml'])
            with open(pointer[1]+'.xml','r') as f:
                result = f.read()
            '''
        bnglFile = libsbml2bngl.readFromString(xmlFile,
                                               reaction, False, None, True)

        result = annotationExtender.expandAnnotation(pointer[1], bnglFile)
        self.addToDict(ticket, result)
        print 'success',

    def generateGraph(self, ticket, bnglContents, graphtype):
        print ticket
        pointer = tempfile.mkstemp(suffix='.bngl', text=True)
        with open(pointer[1], 'w') as f:
            f.write(bnglContents)
        try:
            if graphtype in ['regulatory', 'contactmap']:
                consoleCommands.setBngExecutable(bngDistro)
                consoleCommands.generateGraph(pointer[1], graphtype)
                name = pointer[1].split('.')[0].split('/')[-1]
                with open('{0}_{1}.gml'.format(name, graphtype), 'r') as f:
                    graphContent = f.read()

                gml = networkx.read_gml('{0}_{1}.gml'.format(name, graphtype))
                result = gml2cyjson(gml, graphtype=graphtype)
                jsonStr = json.dumps(result, indent=1, separators=(',', ': '))

                result = {'jsonStr': jsonStr, 'gmlStr': graphContent}
                self.addToDict(ticket, result)
                os.remove('{0}_{1}.gml'.format(name, graphtype))
                print 'success', ticket

            elif graphtype in ['sbgn_er']:
                consoleCommands.setBngExecutable(bngDistro)
                consoleCommands.generateGraph(pointer[1], 'contactmap')
                name = pointer[1].split('.')[0].split('/')[-1]
                # with open('{0}_{1}.gml'.format(name,'contactmap'),'r') as f:
                #   graphContent = f.read()
                graphContent = networkx.read_gml(
                    '{0}_{1}.gml'.format(name, 'contactmap'))
                sbgn = libsbgn.createSBNG_ER_gml(graphContent)
                self.addToDict(ticket, sbgn)
                os.remove('{0}_{1}.gml'.format(name, 'contactmap'))
                print 'success', ticket
            elif graphtype in ['std']:
                consoleCommands.setBngExecutable(bngDistro)
                consoleCommands.bngl2xml(pointer[1])
                xmlFileName = pointer[1].split('.')[0] + '.xml'
                xmlFileName = xmlFileName.split(os.sep)[-1]

                gmlGraph = stdgraph.generateSTDGML(xmlFileName)
                os.remove('{0}.gml'.format(xmlFileName))
                #result = gml2cyjson(gmlGraph, graphtype=graphtype)
                #jsonStr = json.dumps(result, indent=1, separators=(',', ': '))

                #result = {'jsonStr': jsonStr, 'gmlStr': gmlGraph}

                self.addToDict(ticket, ''.join(gmlGraph))
                print 'success', ticket
        except:
            self.addToDict(ticket,-5)
            print 'failure',ticket
        finally:
            task.deferLater(reactor, 600,  freeQueue, ticket)

    def xmlrpc_generateGraph(self, bbnglFile, graphtype):
        counter = next_id()
        bnglFile = bbnglFile.data
        reactor.callInThread(self.generateGraph, counter, bnglFile, graphtype)
        processDict[counter] = -2
        return counter
        # self.generateGraph(counter,bnglFile,graphtype)
        # return counter

    def xmlrpc_generateAnnotations(self, bxmlFile):
        counter = next_id()
        xmlFile = bxmlFile.data
        reactor.callInThread(self.generateAnnotation, counter, xmlFile)
        # reactor.callInThread(self.atomize,counter,xmlFile,False)

        processDict[counter] = -2
        return counter

    def xmlrpc_atomize(self, bxmlFile, atomize=False, reaction='config/reactionDefinitions.json', species=None, buser = None):
        counter = next_id()
        xmlFile = bxmlFile.data
        if buser:
            user = buser.data
        else:
            user = None
        reactor.callInThread(self.atomize, counter, xmlFile, atomize, user)
        processDict[counter] = -2
        # result = threads.deferToThread(libsbml2bngl.readFromString,xmlFile,
        #                                     reaction,True,None,atomize)
        # result = libsbml2bngl.readFromString(xmlFile,
        #                                     reaction,True,None,atomize)

        return counter

    def xmlrpc_getMoleculeTypes(self, bbnglFile, bbnglFile2):
        """receives bbnglFile returns the molecule types in the BNGL"""
        counter = next_id()
        bnglFile = bbnglFile.data
        bnglFile2 = bbnglFile2.data
        reactor.callInThread(self.extractMoleculeTypes, counter, bnglFile, bnglFile2)
        # process is ready to start status
        processDict[counter] = -2

        return counter

    def xmlrpc_compareFiles(self, bbnglFile, bbnglFile2, mappingScript):
        """receives bbnglFile returns the molecule types in the BNGL"""
        counter = next_id()
        bnglFile = bbnglFile.data
        bnglFile2 = bbnglFile2.data
        mappingDict = yaml.load(mappingScript.data)
        print mappingDict, type(mappingDict)
        reactor.callInThread(self.compareFiles, counter, bnglFile, bnglFile2, mappingDict)
        # process is ready to start status
        processDict[counter] = -2

        return counter

    def xmlrpc_getDict(self, ticketNumber):
        if ticketNumber in processDict:
            #return processDict.pop(ticketNumber)
            return processDict[ticketNumber]
        else:
            return -1



    def xmlrpc_isready(self, ticketNumber):
        if ticketNumber in processDict:
            if processDict[ticketNumber] == -2:
                return -2
            return 1
        else:
            return -1


if __name__ == '__main__':
    r = AtomizerServer()
    reactor.listenTCP(9000, server.Site(r))
    reactor.run()
    #f = open('XMLExamples/curated/BIOMD0000000019.xml')
    #s = f.read()
    # r.generateAnnotation(1,s)
