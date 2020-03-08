
''' Chemtools setup script'''

from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys


class Tox(TestCommand):

    # user_options = [('tox-args=', 'a', "Arguments to pass to tox")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.tox_args = None

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import tox
        import shlex
        args = self.tox_args
        if args:
            args = shlex.split(self.tox_args)
        errno = tox.cmdline(args=args)
        sys.exit(errno)


def readme():
    '''Return the contents of the README.rst file.'''
    with open('README.rst', encoding='utf8') as freadme:
        return freadme.read()


setup(
    author="Lukasz Mentel",
    author_email="lmmentel@gmail.com",
    description="Python tools for quantum chemical calculations",
    include_package_data=True,
    install_requires=[
        'sqlalchemy',
        'scipy >= 0.11',
        'numpy >= 1.7',
        'mendeleev',
        'numba',
        'pandas'],
    entry_points={
        'console_scripts': [
            'writeorbinp = chemtools.calculators.gamessus:writeorbinp',
            'bsconvert = chemtools.cli:bsconvert',
            'bsprint = chemtools.cli:bsprint',
        ]
    },
    license='MIT',
    long_description=readme(),
    long_description_content_type='text/x-rst',
    name='chemtools',
    packages=['chemtools', 'chemtools/pescan', 'chemtools/calculators'],
    scripts=[
        'chemtools/submitgamess.py',
        'chemtools/submitmolpro.py'],
    url='https://github.com/lmmentel/chemtools/',
    version='0.9.1',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'],
    keywords=['basis set', 'optimization', 'quantum chemistry',
              'molecular physics'],
    tests_require=['tox', 'pytest'],
    cmdclass={'test': Tox},
)
