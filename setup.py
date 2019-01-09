from setuptools import setup, Extension, find_packages
from glob import glob
import numpy

setup(
    name='pypropagate',
    version='1.2.5',
    description='A python based paraxial wave propagation framework',

    author='Lars Melchior',
    author_email='lars.melchior@gmail.com',

    packages=find_packages(exclude=['tests*']),

    install_requires=[
        'matplotlib',
        'numpy',
        'scipy',
        'expresso[pycas]'
    ],

    zip_safe=False,

    classifiers=[
        'Programming Language :: Python :: 2.7'
    ],

    ext_modules=[
        Extension('_pypropagate',
                  sources = ['source/python.cpp'],
                  include_dirs=['/usr/include/eigen3/'],
                  library_dirs=['/'],
                  extra_compile_args=['-g', '-std=c++14', '-Wall', '-Wpedantic', '-ffast-math', '-O3']
                  ),
        ],

    test_suite='nose.collector',
    tests_require=['nose','scipy'],
)
