from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='biogrinder',
    version='0.1.0',
    author='Phil Charron',
    author_email='phil.charron@inspection.gc.ca',
    description='Biogrinder is a versatile open-source bioinformatic tool to create simulated omic shotgun and amplicon sequence libraries.',
    url='https://github.com/philcharron-cfia/biogrinder',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'biogrinder=main:main',
        ],
    },
    license='MIT',
    keywords='bioinformatics tools sequence-analysis',
    python_requires='>=3.7'
)
