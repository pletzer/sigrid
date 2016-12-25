
import sys
import sigrid
import glob

# set some package attributes
# __version__ = 

# what packages should be imported if the user says from sigrid import *
#__all__ = ['conserveInterp2D']

# open the shared library when importing this module
# on Darwin the suffix is also .so!
from ctypes import CDLL
suffix = 'so'
if sys.platform == 'win32':
    suffix = 'dll'

# try opening the shared library. Depending on the python
# version the library might be called libsigrid.cpython-35m-...so
# keep on trying...
libs = glob.glob(sigrid.__path__[0] + '/*sigrid*.' + suffix)
so = None
for lib in libs:
    try:
        so = CDLL(lib)
        break
    except:
        pass
if not so:
    print('Unable to open shared library libsigrid!')

