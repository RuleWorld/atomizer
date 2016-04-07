import shutil
import platform
import os

version = '1.1'
destdir = os.path.join(version, '{0}-{1}'.format(platform.system(), platform.architecture()[0]))
try:
    os.makedirs(destdir)
except OSError:
    pass
print 'moving Atomizer to {0}\n'.format(destdir)
if platform.system() != 'Windows':
    shutil.move('sbmlTranslator', os.path.join(destdir, 'sbmlTranslator'))
else:
    shutil.move('sbmlTranslator.exe', os.path.join(destdir, 'sbmlTranslator.exe'))