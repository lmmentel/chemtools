#!/usr/bin/env python

from distutils.core import setup, Extension

setup(
    name='chemtools',
    version='0.2.0',
    description = "Python tools for quantum chemical calculations",
    author = "Lukasz Mentel",
    author_email = "lmmentel@gmail.com",
    py_modules = ['basisopt', 'basisset', 'dmft', 'gamessus', 'molecule', 'molpro'],
    license = open('LICENSE.txt').read(),
    long_description = open('README.txt').read(),
    url = 'https://bitbucket.org/lukaszmentel/chemtools/',
)
