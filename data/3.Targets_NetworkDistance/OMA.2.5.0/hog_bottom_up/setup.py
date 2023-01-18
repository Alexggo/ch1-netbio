from setuptools import setup, find_packages
import os
from io import open


name = 'gethogs'
requirements = ['lxml', 'numpy', 'biopython']


__version__ = "Undefined"
for line in open('{}/__init__.py'.format(name.lower())):
    if line.startswith('__version__'):
        exec(line.strip())

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=name,
    version=__version__,
    author='Dessimoz Lab - Laboratory of Computational Evolutionary Biology and Genomics',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='A tool to infer Hierarchical Orthologous Groups (HOGs) from pairwise orthologs',
    long_description=long_description,
    keywords=['orthology', 'HOGs', 'orthoxml', 'algorithm'],
    url='https://github.com/DessimozLab/gethogs',
    license='MIT',
    classifiers=[
         'Development Status :: 4 - Beta',
         'Environment :: Console',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'License :: OSI Approved :: MIT License',
         'Programming Language :: Python :: 2',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3',
         'Programming Language :: Python :: 3.5',
         'Programming Language :: Python :: 3.6',
         'Programming Language :: Python :: 3.7',
         'Programming Language :: Python :: 3.8',
         'Programming Language :: Python :: 3.9',
         ],
    scripts=['bin/warthogs.py'],
    packages=find_packages(),
    install_requires=requirements,
)
