#!/usr/bin/env python

import unittest
from unittest import mock

from maflib.column import MafColumnRecord

from aliquotmaf.converters import collection as MOD


class ThisTestCase(unittest.TestCase):
    pass


class Test_InputData(ThisTestCase):
    CLASS_OBJ = MOD.InputData

    def setUp(self) -> None:
        super().setUp()
        self.mock_value = "foo"
        self.mock_col = "bar"

    def test_obj_creation(self):
        found = self.CLASS_OBJ(self.mock_value, self.mock_col)

        assert self.mock_value in str(found)
        assert self.mock_col in str(found)
        assert found.transformed is None

    @mock.patch("aliquotmaf.converters.collection.get_builder")
    def test_InputData_build(self, builder_mock):
        mock_record = mock.MagicMock(spec_set=MafColumnRecord)
        mock_scheme = "scheme"
        builder_mock.return_value = mock_record

        found = self.CLASS_OBJ(self.mock_value, self.mock_col)

        found.__build__(mock_scheme)


class Test_InputCollection(ThisTestCase):
    def setUp(self) -> None:
        super().setUp()

    def test_object_creation(self):
        obj = MOD.InputCollection()

        self.assertIsInstance(obj, MOD.InputCollection)

    def test_add(self):
        pass


# __END__
