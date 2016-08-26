from distutils.core import setup, Extension 
from Cython.Distutils import build_ext
from Cython.Build import cythonize 

import numpy as np

ext_modules = Extension(name="create_depth_array",
                        sources=["create_depth_array.pyx"],
                        language="c",
                        include_dirs=[np.get_include()]
                        )

setup(
    name = "create_depth_array",
    ext_modules = cythonize(ext_modules),  # accepts a glob pattern
)
