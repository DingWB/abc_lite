from setuptools import setup, find_packages
from setuptools.command.install import install
import glob

script_files = glob.glob('bin/*.py')
setup(
    name='abc_lite',
    version='0.1',
    package_dir={'': 'src'},  # Specifies that packages are under src directory
    packages=find_packages(where='src'),     # Tells setuptools to find packages under src
    scripts=script_files,
    install_requires=[
    ],
)