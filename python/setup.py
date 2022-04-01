#
# Following Jaime's c++ course:
#
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy
from distutils import sysconfig
#import numpy.distutils.intelccompiler
import numpy.distutils.ccompiler
import platform as plt
import sys
import pathlib

#os.system('')
p = pathlib.Path(sys.executable)
root_dir = str(pathlib.Path(*p.parts[0:-2]))


# MacOsx case
if(plt.system() == 'Darwin'):
    root_dir = '/opt/local/'
    CC = 'clang'
    CXX= 'clang++'
    link_opts = ["-bundle","-undefined","dynamic_lookup", "-fopenmp"]

else: # Linux
    #root_dir = '/usr/'
    CC = 'gcc'
    CXX= 'g++'
    link_opts = ["-shared", "-fopenmp"]

os.environ["CC"] = CC
os.environ["CXX"] = CXX

from distutils import sysconfig
sysconfig.get_config_vars()['CC'] =  CC
sysconfig.get_config_vars()['CXX'] = CXX
sysconfig.get_config_vars()['LDSHARED'] = CC
sysconfig.get_config_vars()['CPP'] = CXX


extension_name = "pymirrors"

comp_flags=['-O3','-std=c++14','-march=native','-fPIC','-fopenmp', '-I./src']

extension = Extension(extension_name,
                      sources=["mirror_lib.pyx"], 
                      include_dirs=["./",numpy.get_include(), root_dir+'/include/'],
                      language="c++",
                      extra_compile_args=comp_flags,
                      extra_link_args=comp_flags+link_opts,
                      library_dirs=['./',root_dir+'/lib/'],
                      libraries=[])

extension.cython_directives = {'language_level': "3"}

setup(
    name = extension_name,
    version = '0.0',
    author = 'apy',
    ext_modules=[extension],
    cmdclass = {'build_ext': build_ext}
)

