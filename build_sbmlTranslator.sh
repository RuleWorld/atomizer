#/usr/bin/bash

cd SBMLparser
python3 -m venv atomizer_venv
source atomizer_venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
pip install --target=. python-libsbml
python3 -m PyInstaller sbmlTranslator.spec
deactivate
