# -*- coding: utf-8 -*-
# @Time    : 2024/12/10 15:01
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_bar.py

from hastat.viz.bar import HapBar
import tomli

with open('bar.toml', 'rb') as f:
    config = tomli.load(f)

bar = HapBar(config)
bar.plot()
