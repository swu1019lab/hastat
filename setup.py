# -*- coding: utf-8 -*-
# @Time    : 2024/3/30 10:09
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : setup.py.py

from setuptools import setup, find_packages

setup(
    name='hastat',
    version='0.0.2',
    packages=find_packages(),
    url='https://github.com/swu1019lab/hastat',
    license='BSD License',
    author='Xiaodong Li',
    author_email='lxd1997xy@163.com',
    description='A library for statistical analysis of genotype data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'numpy>=1.21.2',
        'pandas>=1.3.3',
        'scipy>=1.7.1',
        'statsmodels>=0.13.0',
        'scikit-allel>=1.3.6',
        'gffutils>=0.10.1',
        'pysam>=0.17.0',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'hastat=hastat.hastat:main',
        ],
    },
)