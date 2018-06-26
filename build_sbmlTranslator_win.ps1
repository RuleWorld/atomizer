cd .\SBMLparser
# You must add Python to your path for this to work.
# This is pretty easy to do with Anaconda
python -m venv .\atomizer_venv
.\atomizer_venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -r .\requirements_win.txt
pip install --target=. python-libsbml
pyinstaller .\sbmlTranslator.spec
deactivate
cd ..
mkdir bin
cp  .\SBMLparser\dist\sbmlTranslator.exe .\bin
