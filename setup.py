#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    author = "Lukasz Mentel",
    author_email = "lmmentel@gmail.com",
    description = "Python tools for quantum chemical calculations",
    install_requires = [
        'periodic >= 2.0',
        'scipy >= 0.11',
        'numpy >= 1.7',
    ],
    license = open('LICENSE.txt').read(),
    long_description = readme(),
    name = 'chemtools',
    packages = ['chemtools'],
    url = 'https://bitbucket.org/lukaszmentel/chemtools/',
    version = '0.2.4',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords = 'basis set optimization quantum chemistry molecular physics',
)
