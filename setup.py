''' Chemtools '''

from setuptools import setup

def readme():
    '''Return the contents of the README.rst file.'''
    with open('README.rst') as freadme:
        return freadme.read()

#    install_requires = [
#        'scipy >= 0.11',
#        'numpy >= 1.7',
#    ],

setup(
    author = "Lukasz Mentel",
    author_email = "lmmentel@gmail.com",
    description = "Python tools for quantum chemical calculations",
    include_package_data = True,
    license = open('LICENSE.rst').read(),
    long_description = readme(),
    name = 'chemtools',
    packages = ['chemtools'],
    url = 'https://bitbucket.org/lukaszmentel/chemtools/',
    version = '0.3.1',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: MIT',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords = 'basis set optimization quantum chemistry molecular physics',
)


