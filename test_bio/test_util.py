"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

from unittest import TestCase, main

from ruse.bio.bio_util import is_dna, is_protein


class UtilTestCase(TestCase):

    def test_is_dna(self):
        self.assertFalse(is_dna("EVRNAK"))
        self.assertTrue(is_dna("GCTAA"))

    def test_is_protein(self):
        self.assertTrue(is_protein("EVRNAK"))
        self.assertTrue(is_protein("GCTAA"))
        self.assertTrue(is_protein("ABCDEGGHIJKLMNOPQRSTUVWXYZ-*"))
        self.assertFalse(is_protein("123456"))


if __name__ == "__main__":
    main()
