''' Chemtools setup script'''

from setuptools.command.test import test as TestCommand
import sys

class Tox(TestCommand):

    #user_options = [('tox-args=', 'a', "Arguments to pass to tox")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.tox_args = None
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import tox
        import shlex
        args = self.tox_args
        if args:
            args = shlex.split(self.tox_args)
        errno = tox.cmdline(args=args)
        sys.exit(errno)


from setuptools import setup

def readme():
    '''Return the contents of the README.rst file.'''
    with open('README.rst') as freadme:
        return freadme.read()

setup(
    author = "Lukasz Mentel",
    author_email = "lmmentel@gmail.com",
    description = "Python tools for quantum chemical calculations",
    include_package_data = True,
    install_requires = [
        'sqlalchemy',
        'scipy >= 0.11',
        'numpy >= 1.7',
        'mendeleev',
    ],

    license = open('LICENSE.rst').read(),
    long_description = readme(),
    name = 'chemtools',
    packages = ['chemtools', 'chemtools/pescan'],
    scripts = [
        'chemtools/submitgamess.py',
        'chemtools/submitmolpro.py',
        'scripts/bsprint',
        'scripts/bsconvert',
               ],
    url = 'https://bitbucket.org/lukaszmentel/chemtools/',
    version = '0.5.0',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: MIT',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords = ['basis set', 'optimization', 'quantum chemistry', 'molecular physics'],
    tests_require=['tox', 'pytest'],
    cmdclass = {'test': Tox},
)
