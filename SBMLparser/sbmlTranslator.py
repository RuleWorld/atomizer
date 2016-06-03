# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:38:21 2013

@author: proto
"""

import libsbml2bngl as ls2b
import argparse
import yaml

def defineConsole():
    parser = argparse.ArgumentParser(description='SBML to BNGL translator')
    parser.add_argument('-i', '--input-file', type=str, help='input SBML file', required=True)
    parser.add_argument('-t', '--annotation', action='store_true', help='keep annotation information')
    parser.add_argument('-o', '--output-file', type=str, help='output SBML file')
    parser.add_argument('-c', '--convention-file', type=str, help='Conventions file')
    parser.add_argument('-n', '--naming-conventions', type=str, help='Naming conventions file')
    parser.add_argument('-u', '--user-structures', type=str, help='User defined species')
    parser.add_argument('-id', '--molecule-id', action='store_true', help='use SBML molecule ids instead of names. IDs are less descriptive but more bngl friendly. Use only if the generated BNGL has syntactic errors')
    parser.add_argument('-nc','--no-conversion', action='store_true', help='do not convert units. Copy straight from sbml to bngl')
    parser.add_argument('-a', '--atomize', action='store_true', help='Infer molecular structure')
    parser.add_argument('-p', '--pathwaycommons', action='store_true', help='Use pathway commons to infer molecule binding. This setting requires an internet connection and will query the pathway commons web service.')
    parser.add_argument('-b', '--bionetgen-analysis', type=str, help='Set the BioNetGen path for context post analysis.')
    parser.add_argument('-s','--isomorphism-check', action='store_true', help='disallow atomizations that produce the same graph structure')
    parser.add_argument('-I','--ignore', action='store_true', help='ignore atomization translation errors')

    return parser


def checkInput(namespace):
    options = {}
    options['inputFile'] = namespace.input_file

    conv, useID, naming = ls2b.selectReactionDefinitions(options['inputFile'])
    options['outputFile'] = namespace.output_file if namespace.output_file is not None else options['inputFile'] + '.bngl'
    options['conventionFile'] = namespace.convention_file if namespace.convention_file is not None else conv
    options['userStructure'] = namespace.user_structures
    options['namingConventions'] = namespace.naming_conventions if namespace.naming_conventions is not None else naming
    options['useId'] = namespace.molecule_id
    options['annotation'] = namespace.annotation
    options['atomize'] = namespace.atomize
    options['pathwaycommons'] = namespace.pathwaycommons
    options['bionetgenAnalysis'] = namespace.bionetgen_analysis
    options['isomorphismCheck'] = namespace.isomorphism_check
    options['ignore'] = namespace.ignore
    options['noConversion'] = namespace.no_conversion
    return options


def main():
    parser = defineConsole()
    namespace = parser.parse_args()

    options = checkInput(namespace)
    returnArray = ls2b.analyzeFile(options['inputFile'], options['conventionFile'], options['useId'], options['namingConventions'],
                                   options['outputFile'], speciesEquivalence=options['userStructure'],
                                   atomize=options['atomize'], bioGrid=False, pathwaycommons=options['pathwaycommons'], ignore=options['ignore'], noConversion = options['noConversion'])

    if namespace.bionetgen_analysis and returnArray:
        ls2b.postAnalyzeFile(options['outputFile'], namespace.bionetgen_analysis, returnArray.database)

    if namespace.annotation and returnArray:
        with open(options['outputFile'] + '.yml', 'w') as f:
            f.write(yaml.dump(returnArray.annotation, default_flow_style=False))

if __name__ == "__main__":
    main()
