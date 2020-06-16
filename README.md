Project for the translation of SBML files into BNGL files. You can find more
information about the project [here](https://ruleworld.github.io/atomizer/).

[![Build Status](https://travis-ci.org/RuleWorld/atomizer.svg?branch=master)](https://travis-ci.org/RuleWorld/atomizer) [![Build status](https://ci.appveyor.com/api/projects/status/bst5csm5vyo1q0vp?svg=true)](https://ci.appveyor.com/project/jczech/atomizer)

## Download Latest Test Builds

These builds are the from the head of master and are not guaranteed to be
stable. Use at your own risk.

* [Linux](https://bintray.com/jczech/atomizer/download_file?file_path=sbmlTranslator_linux)
* [OSX](https://bintray.com/jczech/atomizer/download_file?file_path=sbmlTranslator_osx)
* [Windows x64](https://bintray.com/jczech/atomizer/download_file?file_path=sbmlTranslator_win.zip)

## Requirements

You should use Python 3.5 or later. A number of Python libraries need to be
installed, but these can be automatically installed in a Python virtual
environment, which is exactly what the makefile does under Linux and OSX or
what the PowerShell script does for Windows.

### Ubuntu

    sudo apt-get install -y python3-dev
    sudo apt-get install -y python3-venv

### Windows 10

You'll need to install Python 3. We recommend using Anaconda for this.

### Optional for Developers

If you want to run `./atomizer/SBMLparser/sbmlTranslator.py` directly, you can
install all the Python requirements at the system level by doing:

    cd SBMLparser
    sudo pip3 install -r requirements.txt

Alternatively, you can do this at the user level:

    cd SBMLparser
    pip3 install --user -r requirements.txt

## Installation 

### OSX and Linux

From the top-level directory, type the following:

    make; make install

This will create an `sbmlTranslator` binary in `./atomizer/bin`. This is made
using PyInstaller.

### Installation for Windows

Right click on `./atomizer/build_sbmlTranslator_win.ps1` and select `Run with
PowerShell`.

## Execution

    ./sbmlTranslator -i /path/to/sbml/file [-a] [-o output.bngl]

Optional arguments
- [-a] activates the atomizer. Otherwise only a flat translation will be provided (no molecular structure)
- [-p] activates pathway commons querying. An internet connection is required
- [-b <path/to/bionetgen>] Enables post atomization analysis.
- See [-h] for a full list of arguments


Directory Structure:
- SBMLparser: Directory containing the main project.
    - sbmlTranslator.py: Entry level file
    - atomizer: Contains the main implicit assumption extraction code
        - analyzeSBML.py: lexical analysis engine
        - detectOntology.py: helper file for the lexical analysis engine
          (levenshtein distance, convention discovery, etc.)
        - analyzeRDF.py: extracts annotation information and gets it into the
          atomizer data model
        - resolveSCT.py: Creates the species composition table using lexical
          analysis, annotation information and the stoichiometry matrix
        - moleculeCreation.py: Creates the structured molecules for use in
          graph-creation

    - utils: contains several utility scripts used during the atomization process
        - consoleCommands.py: Interfaces with the BioNetGen console
        - pathwaysCommons.py: Queries BioGrid and pathways commons given an RDF annotation entry
- stats: Extracts statistics of a set of atomized BioModels files
