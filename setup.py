import sys
import glob

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

sources = glob.glob('src/*.cpp')
try:
    sources.remove('src/main.cpp')
except ValueError:
    pass

setup(name='Python wrapper around SIMD partial order alignment library',
    ext_modules=cythonize([Extension(
        'spoapy',
        sources=['spoapy/spoapy.pyx'] + sources,
        include_dirs=['include/', 'include/spoa', 'src'],
        language='c++',
        extra_compile_args=['-std=c++11'],
        extra_link_args=['-std=c++11']
    )],
    compiler_directives={'language_level' : sys.version_info[0]}))
