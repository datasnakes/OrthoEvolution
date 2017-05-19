from unittest import TestCase
from Orthologs import Align

class Test(TestCase):
    def test_is_string(self):
        s = Align.joke()
        self.assertTrue(isinstance(s))

help(Align)