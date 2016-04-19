Project for the translation of SBML files into BNGL files 

[![Build Status](https://travis-ci.org/RuleWorld/atomizer.svg?branch=master)](https://travis-ci.org/RuleWorld/atomizer) [![Build status](https://ci.appveyor.com/api/projects/status/rb4sci41f2fy62il?svg=true)](https://ci.appveyor.com/project/jjtapia/atomizer)


## Requirements
libsbml, networkx (for state transition diagram creation), pexpect (for post atomization analysis), 

## Installation:

> make; make install

## Execution

> ./sbmlTranslator -i /path/to/sbml/file [-a] [-o output.bngl]

Optional arguments
- [-a] activates the atomizer. Otherwise only a flat translation will be provided (no molecular structure)
- [-p] activates pathway commons querying. An internet connection is required
- [-b <path/to/bionetgen>] Enables post atomization analysis.
- See [-h] for a full list of arguments


Directory Structure:
|
- SBMLparser: Directory containing the main project.
  | 
    - sbmlTranslator.py: Entry level file
    - atomizer: Contains the main implicit assumption extraction code
    	- analyzeSBML.py: lexical analysis engine
        - detectOntology.py: helper file for the lexical analysis engine (levenshtein distance, convention
        					 discovery, etc.)
        - analyzeRDF.py: extracts annotation information and gets it into the atomizer data model
        - resolveSCT.py: Creates the species composition table using lexical analysis, annotation information and
        				 the stoichiometry matrix
        - moleculeCreation.py: Creates the structured molecules for use in graph-creation
        
    - utils: contains several utility scripts used during the atomization process
    	- consoleCommands.py: Interfaces with the BioNetGen console
        - pathwaysCommons.py: Queries BioGrid and pathways commons given an RDF annotation entry
- stats: Extracts statistics of a set of atomized BioModels files

