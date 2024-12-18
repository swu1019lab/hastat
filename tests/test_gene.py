# -*- coding: utf-8 -*-
# @Time    : 2024/12/10 18:25
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_gene.py

from hastat.viz.gene import HapGene
import tomli

with open('gene.toml', 'rb') as f:
    config = tomli.load(f)

gene = HapGene(config)
gene.plot()
