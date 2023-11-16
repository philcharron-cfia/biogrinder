from setuptools import setup, find_packages

setup(
    name='biogrinder',
    version='0.1.0',
    author='Phil Charron',
    author_email='phil.charron@inspection.gc.ca',
    description='Grinder is a versatile open-source bioinformatic tool to create simulated omic shotgun and amplicon sequence libraries.'
    url='https://github.com/philcharron-cfia/biogrinder'
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [
            'biogrinder = src.main:main',
        ],
    },
    # Add other parameters like author, license, description, etc.
)
