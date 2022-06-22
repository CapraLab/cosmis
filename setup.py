# setup.py

from setuptools import setup

setup(
    name='cosmis',
    version='1.0',
    description='COSMIS is a framework for quantifying the mutational constraint on amino acid sites in 3D spatial neighborhoods.'
    author='Bian Li',
    author_email='bian.li@vanderbilt.edu',
    packages=['cosmis'],
    scripts=[
        'scripts/cosmis',
        'scripts/cosmis-sp'
    ],
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'pytest',
        'wget',
        'zlib'
    ]
)
