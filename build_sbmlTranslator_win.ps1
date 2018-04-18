cd .\SBMLparser
# You must add Python to your path for this to work.
# This is pretty easy to do with Anaconda
python -m venv .\atomizer_venv
.\atomizer_venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -r .\requirements.txt
pip install --target=. python-libsbml
pyinstaller .\sbmlTranslator.spec
deactivate