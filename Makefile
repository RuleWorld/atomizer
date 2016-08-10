

SHELL = /bin/bash
BUILD = ./build
DIST = ./dist

.PHONY: all install clean test


all:
	if ! test -d pyinstaller2 ; then \
		unzip utils/pyinstaller2.zip;   \
	fi ;


ifeq ($(OS),Windows_NT)
ifeq ($(shell uname -o),Cygwin)
	python pyinstaller2/pyinstaller.py utils/sbmlTranslator.spec ;
else
	python pyinstaller2/pyinstaller.py utils/sbmlTranslator_windows.spec ;
endif
else
	virtualenv --no-site-packages venv
	source venv/bin/activate
	pip install --user -r requirements.txt
	python pyinstaller2/pyinstaller.py utils/sbmlTranslator.spec ;
	deactivate
endif

install:
	mkdir -p bin
ifeq ($(OS),Windows_NT)
    ifeq ($(shell uname -o),Cygwin)
	    cp  ${DIST}/sbmlTranslator bin/sbmlTranslator.exe;
    else
	    cp  ${DIST}/sbmlTranslator.exe bin/sbmlTranslator.exe;
    endif
else
	cp  ${DIST}/sbmlTranslator bin/sbmlTranslator;
endif


test:
	cd SBMLparser; PYTHONPATH=$PYTHONPATH:. nosetests --with-doctest; cd ..;
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
	if test -d venv; then \
		rm -rf venv;	\
	fi

	find . -name '*.pyc' -delete

	
