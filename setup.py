#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: domenico.somma@glasgow.ac.uk
"""

from setuptools import setup

setup(
    name='papillon',
    version='0.1.1',
    py_modules=['papillon'],
    description='A Python module to read and plot (cuffdiff) Galaxy RNA-seq data',
    author='Domenico Somma',
    author_email='domenico.somma@glasgow.ac.uk',
    license='Mozilla Public License 2.0',
    url='https://github.com/domenico-somma/Papillon/',
    python_requires='>=3.3, <4',
    install_requires=[
        "pandas >= 0.17.1",
        "Seaborn >= 0.8.1",
    ],
)
