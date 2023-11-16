from setuptools import setup, find_packages
import unittest

def test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('.', pattern='test*.py')
    #test_runner = unittest.TextTestRunner(verbosity=2)
    #test_runner.run(test_suite)
    return test_suite

setup(
    name='biogrinder',
    version='0.1.0',
    author='Phil Charron',
    author_email='phil.charron@inspection.gc.ca',
    description='Grinder is a versatile open-source bioinformatic tool to create simulated omic shotgun and amplicon sequence libraries.',
    url='https://github.com/philcharron-cfia/biogrinder',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[],
    entry_points={
        'console_scripts': [
            'biogrinder = src.main:main',
        ],
    },
    license='MIT',
    keywords='bioinformatics tools sequence-analysis',
    python_requires='>=3.6',
    test_suite='setup.test_suite'
)
