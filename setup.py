# -*- coding: utf-8 -*-
# @Time    : 2024/3/30 10:09
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : setup.py.py

from setuptools import setup, find_packages

setup(
    name='hastat',
    version='0.0.5',
    packages=find_packages(),
    url='https://github.com/swu1019lab/hastat',
    license='BSD License',
    author='Xiaodong Li',
    author_email='lxd1997xy@163.com',
    description='A library for gene haplotypes analysis and statistics.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
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