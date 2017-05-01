# -*- coding: utf-8 -*-
from logit import Logit

foo = Logit("cray", "foo")
bar = Logit("cray", "bar")

foo.scriptinfo()
foo.log.info('hey')
bar.log.info('hey')

