
import sys
import sigrid
import glob

# set some package attributes
# __version__ = 

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
libsigrid = None
for lib in libs:
    try:
        libsigrid = CDLL(lib)
        break
    except:
        pass
if not libsigrid:
    print('Unable to open shared library libsigrid!')

