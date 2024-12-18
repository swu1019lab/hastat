# -*- coding: utf-8 -*-
# @Time    : 2024/12/10 15:51
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_pie.py


from hastat.viz.pie import HapPie
import tomli

with open('pie.toml', 'rb') as f:
    config = tomli.load(f)

pie = HapPie(config)
pie.plot()
