# -*- coding: utf-8 -*-

from glob import glob
from os.path import basename, splitext
from setuptools import find_packages, setup

setup(
    name='vptools',
    version='0.1',
    packages=find_packages(where='vptools'),
    package_dir={'': 'vptools'},
    py_modules=[splitext(basename(path))[0] for path in glob('vptools/*/*.py')],
)