# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 00:46:05 2014

@author: proto
"""

import unittest
import os
import sys
import numpy as np
import pexpect
from subprocess import call

pathname = os.path.abspath(os.path.dirname(__file__))
print '---', pathname
sys.path.insert(0, '.')

sys.path.insert(0, os.path.join(pathname, '..'))
sys.path.insert(0, os.path.join(pathname, '..', 'SBMLparser'))


import SBMLparser.libsbml2bngl as libsbml2bngl

bionetgenBinary = os.path.join(pathname, 'BioNetGen-2.2.6-stable', 'BNG2.pl')
def bnglExecution(bnglFile, settings):
    os.chdir(os.path.join(pathname, 'tmp'))
    try:
        bngconsole = pexpect.spawn('perl {0} --console'.format(bionetgenBinary))
        bngconsole.expect('BNG>')
        bngconsole.sendline('load {0}.bngl'.format(bnglFile))

        bngconsole.expect('BNG>')
        bngconsole.sendline('action generate_network()')
        bngconsole.expect('BNG>')
        bngconsole.sendline('action simulate({{method=>"ode",t_start=>{0},t_end=>{1},n_steps=>{4},atol=>{2},rtol=>{3}}})'.format(
            settings['start'][0], settings['duration'][0], settings['absolute'][0], settings['relative'][0], settings['steps'][0]))

        bngconsole.expect('BNG>')
        bngconsole.close()
        with open('{0}.gdat'.format(bnglFile)) as f:
            output = f.readlines()
        output = [x.strip().split(' ') for x in output]
        header = [x for x in output[0] if x not in ['#', '']]
        values = np.array([[float(y) for y in x if y != '']
                           for x in output[1:]])
        headerIndex = []
        headerName = []
        for x in settings['variables']:
            for idx, ele in enumerate(header):
                if ele.startswith('{0}_'.format(x)):
                    headerIndex.append(idx)
                    headerName.append('_'.join(ele.split('_')[:-1]))
                    break
        #headerIndex = [header.index(x) for x in settings['variables']]
        values = values[:, headerIndex]
    finally:
        os.chdir('..')

    return values, settings['absolute'][0], headerName


def parseCSV(fileName, headers):
    with open(fileName, 'r') as f:
        inputLines = f.readlines()
    inputLines = [x.strip().split(',') for x in inputLines]
    content = np.array([[float(y) for y in x] for x in inputLines[1:]])
    header = [x.strip() for x in inputLines[0]]
    headerIndex = [header.index(x) for x in headers]
    values = content[:, headerIndex]
    return values


class ParametrizedTestCase(unittest.TestCase):

    """ TestCase classes that want to be parametrized should
        inherit from this class.
    """

    def __init__(self, methodName='runTest', param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_klass, param=None):
        """ Create a suite containing all tests taken from the given
            subclass, passing them the parameter 'param'.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, param=param))
        return suite


class AtomizationTestCase(ParametrizedTestCase):

    '''
    Test for ability to ruleify
    '''

    def setUp(self):
        #reactionDefinitions, useID = libsbml2bngl.selectReactionDefinitions('BIOMD%010i.xml' % self.param)
        #spEquivalence = detectCustomDefinitions(bioNumber)
        pass

    def extractSimulationSettings(self, fileName):
        with open(fileName, 'r') as f:
            settings = f.readlines()

        settings = [x.strip().split(':') for x in settings]
        settingsDict = {x[0]: x[1].split(',') for x in settings if len(x) == 2}
        for element in settingsDict:
            settingsDict[element] = [x.strip() for x in settingsDict[element]]
        return settingsDict

    def test_valid(self):

        if (self.param) is None:
            return
        print '+++', self.param
        # deepreload.reload(libsbml2bngl)
        reactionDefinitions = os.path.join(
            pathname, '..', 'reactionDefinitions', 'reactionDefinition7.json')
        namingConventions = os.path.join(
            pathname, '..', 'config', 'namingConventions.json')
        outputFile = os.path.join(
            pathname, 'tmp', 'output{0}.bngl'.format(self.param[1]))
        #libsbml2bngl.analyzeFile('{0}/{1}/{1}-sbml-l2v4.xml'.format(self.param[0], self.param[1]), reactionDefinitions,
        #                         False, namingConventions,
        #                         outputFile=outputFile, speciesEquivalence=None, atomize=True, bioGrid=False)
        call([os.path.join(pathname,'..','dist','sbmlTranslator'),'-i',os.path.join(pathname,self.param[0],self.param[1],'{0}-sbml-l2v4.xml'.format(self.param[1])),
                           '-o',outputFile,'-a','-nc'])
        settings = self.extractSimulationSettings(os.path.join(pathname,self.param[0],self.param[1],'{0}-settings.txt'.format(self.param[1])))

        bnglValues, atol, validHeaders = bnglExecution('output{0}'.format(self.param[1]), settings)
        referenceValues = parseCSV(
            '{0}/{1}/{1}-results.csv'.format(self.param[0], self.param[1]), validHeaders)
        print '---', float(atol), (((bnglValues - referenceValues)**2).mean())**0.5
        self.assertAlmostEqual(
            (((bnglValues - referenceValues)**2).mean())**0.5, 0, delta=float(atol))
        dirs = [f for f in os.listdir(
            os.path.join(pathname, 'tmp')) if self.param[1] in f]
        for element in dirs:
            os.remove(os.path.join(pathname, 'tmp', element))

    def tearDown(self):
        dirs = [f for f in os.listdir(
            os.path.join(pathname, 'tmp')) if not f.endswith('bngl')]
        for element in dirs:
            os.remove(os.path.join(pathname, 'tmp', element))

AtomizationTestCase.slow = 1

class TestValid(ParametrizedTestCase):

    '''
    Test for whether a file is recognized as correct by bng --check
    '''
    # skip stoichiometry math
    xdirs = ['00068', '00069', '00070', '00518']
    # non supported operands
    xdirs.extend(['00028', '00201', '00269', '00279', '00949', '00951'])
    xdirs.append('00588')
    xdirs.extend(['00922','00949','00951'])
    #weird lambda functions
    xdirs.extend(['00834', '01034'])
    # time functions not fully supported
    xdirs.extend(['00853', '00856', '00859', '00862', '00896'])
    dirs = [f for f in os.listdir(os.path.join(pathname, 'semantic'))]
    dirs.sort()
    dirs = [x for x in dirs if x not in xdirs]
    #dirs=['00465']
    #dirs = ['00813', '00834', '00853', '00856', '00859', '00862', '00896', '01034', '01059']
    #dirs = ['01059']
    #dirs = ['00076', '00077', '00603', '00602']
    #dirs = ['00001']
    suite = unittest.TestSuite()
    for t in [x for x in dirs if x not in xdirs]:
        suite.addTest(ParametrizedTestCase.parametrize(
            AtomizationTestCase, param=[os.path.join(pathname, 'semantic'), t]))

    result =  unittest.TextTestRunner(verbosity=2).run(suite)
    ret = not (result.failures is [] or result.errors is [])
    ret = 0 if ret else 1
    print ret
    sys.exit(ret)    
