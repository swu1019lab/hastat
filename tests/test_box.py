# -*- coding: utf-8 -*-
# @Time    : 2024/12/10 16:34
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_box.py

from hastat.viz.box import HapBox
import tomli

with open('box.toml', 'rb') as f:
    config = tomli.load(f)

box = HapBox(config)
box.plot()
