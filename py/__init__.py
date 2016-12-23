
import sys
import sigrid
import glob

# __version__ = 

# open the shared library when importing this module
# on Darwin the suffix is also .so!
from ctypes import CDLL
suffix = 'so'
if sys.platform == 'win32':
    suffix = 'dll'

# try opening the shared library. Depending on the python
# version the library might be called libcf.cpython-35m-...so
# keep on trying...
libs = glob.glob(sigrid.__path__[0] + '/*sigrid*.' + suffix)
success = False
for lib in libs:
    try:
        print('trying to open {}'.format(lib))
        lib = CDLL(lib)
        success = True
    except:
        pass
    if success:
        break
if not success:
    print('Unable to open shared library!')

