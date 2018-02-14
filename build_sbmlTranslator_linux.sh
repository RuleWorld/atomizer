#/usr/bin/bash

virtualenv2 --no-site-packages venv
source venv/bin/activate
pip2 install -r requirements.txt
python2 pyinstaller2/pyinstaller.py utils/sbmlTranslator.spec
deactivate
