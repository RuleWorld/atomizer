#/usr/bin/bash

#virtualenv2 --no-site-packages venv
cd SBMLparser
python -m venv atomizer_venv
source atomizer_venv/bin/activate
pip install -r requirements.txt
pip install --target=. python-libsbml
pyinstaller sbmlTranslator.spec
deactivate
