SHELL = /bin/bash
BUILD = ./SBMLparser/build
DIST = ./SBMLparser/dist

.PHONY: all install clean test


all:
	./build_sbmlTranslator.sh


install:
	mkdir -p bin
	cp  ${DIST}/sbmlTranslator bin/sbmlTranslator;

test:
	cd SBMLparser; PYTHONPATH=$(PYTHONPATH):. nosetests --with-doctest; cd ..;
ifeq ($(shell uname -s),Linux)
	cd test; tar xfj bionetgen-2.2.6.tar.bz2; tar xfj testsuite.tar.bz2; python testSuite.py;
endif

clean:
	if test -d ${BUILD} ; then \
	    rm -rf ${BUILD} ;          \
	fi ;
	if test -d ${DIST} ; then \
	    rm -rf ${DIST} ;          \
	fi ;
	if test -d pyinstaller2 ; then \
	    rm -rf pyinstaller2 ;          \
	fi ;	
	if test -d test/BioNetGen-2.2.6-stable; then \
		rm -rf test/BioNetGen-2.2.6-stable;	\
	fi
	if test -d SBMLparser/atomizer_venv; then \
		rm -rf SBMLparser/atomizer_venv;	\
	fi
	if test -d SBMLparser/libsbml; then \
		rm -rf SBMLparser/libsbml;	\
	fi
	if test -d SBMLparser/python_libsbml-*.dist-info; then \
		rm -rf SBMLparser/python_libsbml-*.dist-info;	\
	fi

	find . -name '*.pyc' -delete
	find . -name '__pycache__' -delete
