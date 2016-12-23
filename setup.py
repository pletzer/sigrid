#
# Python installation script
#
# Python installation script. Will install the shared library 
# in the directory where python puts its modules (generally 
# /usr/lib/python<VERSION>/site-packages. 
#

from __future__ import print_function
from setuptools import setup, Extension
import glob
import os.path
import sys
import subprocess
import re

def getVersion():
  version_major = "NOT-FOUND"
  version_minor = "NOT-FOUND"
  version_patch = "NOT-FOUND"
  for line in open('CMakeLists.txt').readlines():
  	m = re.search(r'VERSION_MAJOR\s+\"?(\d+)\"?', line)
  	if m:
  		version_major = m.group(1)
  	else:
  		m2 = re.search(r'VERSION_MINOR\s+\"?(\d+)\"?', line)
  		if m2:
  			version_minor = m2.group(1)
  		else:
  			m3 = re.search(r'VERSION_PATCH\s+\"?(\d+)\"?', line)
  			if m3:
  				version_patch = m3.group(1)
  return version_major + '.' + version_minor + '.' + version_patch

def getValuesFromOption(opt, cmd):
  """
  Get the values from a given option in a command
  @param opt e.g. -L or -l
  @param cmd e.g. "-L/usr/local/myLib -lmyLib"
  @return list of values (e.g. ['/usr/local/myLib'])
  """
  values = []
  pat = opt + r'\s*([\w\_\/\\]+)'
  m = re.search(pat, cmd)
  while m:
    # found
    val = m.group(1)
    values.append(val)
    # remove the option-value pair
    cmd = re.sub(opt + '\s*' + val, '', cmd)
    m = re.search(pat, cmd)
  return values

def parseLinkCommand(cmd):
  """
  Find all the library paths and library names from the link command
  @param cmd e.g. "-L/usr/local/myLib -lmyLib"
  @return list of library directories and list of library names
  """
  libdirs = getValuesFromOption('-L', cmd)
  libs = getValuesFromOption('-l', cmd)
  return libdirs, libs

def breakLibraryPath(path):
  """
  Break a library path into -L and -l parts
  @param path e.g. /usr/lib/liblapack.so
  @return library directory (e.g. '/usr/lib') and library name (e.g. 'lapack')
  """
  dirname, libname = os.path.split(path)
  # remove the suffix
  libname = re.sub(r'\..*$', '', libname)
  # remove leading 'lib' if on UNIX
  if sys.platform is not 'win32':
    libname = re.sub(r'^\s*lib', '', libname)
  return dirname, libname


def findLibrary(dirs, name):
  """
  Find library by searching a few common directories
  @param dirs list of directories
  @param name name of the library, without lib and without suffix
  @return full path to the library or '' if not found
  """
  for directory in dirs:
    libname = directory + '/lib' + name + '.*'
    libs = glob.glob(libname)
    if len(libs) > 0:
      # use the first occurrence
      return libs[0]
  return ''

# Get LAPACK and BLAS. If not set then search common locations
dirs = ('/usr/local/lib', '/usr/lib', '/usr/lib64')
lapack_libraries = os.environ.get('LAPACK_LIBRARIES', findLibrary(dirs, 'lapack'
))
blas_libraries = os.environ.get('BLAS_LIBRARIES', findLibrary(dirs, 'blas'))

if not lapack_libraries:
  print('ERROR: could not find lapack -- set environment variable LAPACK_LIBRARIES and rerun')
  sys.exit(1)

if not blas_libraries:
  print('ERROR: could not find blas -- set environment variable BLAS_LIBRARIES and rerun')
  sys.exit(2)

lapack_dir, lapack_lib = breakLibraryPath(lapack_libraries)
blas_dir, blas_lib = breakLibraryPath(blas_libraries)

libdirs = [lapack_dir, blas_dir]
libs = [lapack_lib, blas_lib]

# list all the directories that contain source file to be compiled 
# into a shared library
dirs = ['./cpp',]

# list all the include directories
incdirs = ['./'] + dirs

srcs = []
for d in dirs:
  srcFiles = glob.glob(d + '/*.cpp')
  for s in srcFiles:
    srcs.append(s)

ext_modules = [Extension("libsigrid", # name of the shared library
                          srcs,
                          define_macros=[('HAVE_LAPACK_UNDERSCORE', 1),],
                          include_dirs=incdirs,
                          libraries=libs,
                          library_dirs=libdirs)]

setup(name = "sigrid",
      version = getVersion(),
      long_description = "sgrid -- a package for interpolating fields on structured grids",
      author_email      = "alexander.pletzer@nesi.org.nz",
      url               = "http://github.com/pletzer/sigrid",
      platforms         = ["any"],
      description = "Python access to the sgrid library",
      packages=['sigrid'],
      package_dir={'sigrid': 'py'}, # location of the python files
      ext_modules=ext_modules)
