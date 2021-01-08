#!/usr/bin/env python3

import unittest
from unittest import mock

from aliquotmaf.vcf_to_aliquot.converters import collection as MOD


class ThisTestCase(unittest.TestCase):
    def setUp(self) -> None:
        super().setUp()

    def tearDown(self) -> None:
        super().tearDown()


class Test_InputCollection(ThisTestCase):
    CLASS_OBJ = MOD.InputCollection

    def setUp(self) -> None:
        super().setUp()

    def test_init(self):
        obj = self.CLASS_OBJ()
        self.assertIsInstance(obj, self.CLASS_OBJ)


class Test_InputData(ThisTestCase):
    CLASS_OBJ = MOD.InputData

    def setUp(self) -> None:
        super().setUp()

    def test_init(self):
        value = 'foo'
        column = 'bar'
        obj = self.CLASS_OBJ(value, column)
        self.assertIsInstance(obj, self.CLASS_OBJ)
        self.assertEqual(value, obj.value)
        self.assertEqual(column, obj.column)


# __END__
