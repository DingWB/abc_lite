from setuptools import setup, find_packages
from setuptools.command.install import install


setup(
    name='abc_lite',
    version='0.1',
    package_dir={'': 'src'},  # Specifies that packages are under src directory
    packages=find_packages(where='src'),     # Tells setuptools to find packages under src
    scripts=['bin/run_abc_lite.py'],
    install_requires=[
    ],
)