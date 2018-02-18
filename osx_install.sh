#!/bin/bash
if [ "${TRAVIS_OS_NAME}" = "linux" ]; then
    echo "This should say linux"
    echo $TRAVIS_OS_NAME
    uname -a
    
else
    echo  "This should say osx"
    echo $TRAVIS_OS_NAME
    uname -a
    
    brew update
    brew install python
    brew install libxml2
    brew install libxslt
    brew install libzip
    #brew install lib32z1 lib32z1-dev
    #brew link --overwrite python
    brew unlink python && brew link python
    pip install numpy
    pip install python-libsbml
fi
