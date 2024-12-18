# -*- coding: utf-8 -*-
# @Time    : 2024/12/10 14:46
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : manage.py


class Manager:
    def __init__(self):
        self._registry = {}

    def register(self, name, cls):
        self._registry[name] = cls

    def use(self, name, *args, **kwargs):
        if name not in self._registry:
            raise ValueError(f"The plot type {name} is not supported")
        return self._registry[name](*args, **kwargs)
