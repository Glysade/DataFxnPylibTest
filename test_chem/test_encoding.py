"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

from unittest import TestCase, main

from ruse.chem.chem_data_table_helper import _decode_mol_cell


class TestEncoding(TestCase):

    def test_encoding(self):
        base64_str = 'H4sIAAAAAAAEAKVXTW8TMRC976+YA0jJoa5n/I2gUpXyJRSCWuBA1UNJlxIpDSgNB/4wv4Px2lHT9aI2Jl1Vm+fxy3szY6/3+At6NA3Awezl5Ht7cyA1anRopKWTpgFlQHkACfGD3c3dFUKAzySlbOIgCUnBdHGCpMd4JwWPSphAb+bg1bGgQIt5rkNStSzSaNYNB1KQp0otPMNplVi0t3UsEIRxhtJc5UnXaiEjdcoukt1lef94LV4YjzLN1Z7LV+XIiRBkzq5Ct6tlto8WdE4mR1pqW+fICfTkE4vBf1SaHmKxwgSds+Gkqqy0YRZj0lyPWM2CXue8OLq3ju7yoh5mcZ46FhJGWVmnxTIL6sSiyQ+zPJhdJ5QjF++UIKWr+0VZ6xML0jbP+zrSnN2Q+8UrY+tYFPeadYlFBrW7jvboXRJWUt4xjaRKLcT7JKpco5BU1TiyPmQWJ+1uds8ez4Ii2CBzjXK1ahyhS33P1dJ619Ee/cKVQWcSnzc+1LI4ZWTmk9t9fF9HvB9Yt600UmW/5Bn/rWXLUtzVOIrPxuAqV+Nddgste+VlqwX7WvaodJyBmSVUOwrC2q5jmUUFf+/8snwsC59fXDAhZTcYZWpYsNuZsYcyZAdR1dvImw7SZSxDpow1g7yMuhJlyJcMjIYiNh5WEmrvo7FSZWx85vZR5KNpiRLEPbOPKsDSMUNDv2YAbanMApaOI68vUQ80oMwDhiI7DMWneS+WIcIilqEB3oiWjhmKJ9I+A6Ol4/hrZY0ZUgMaLFCZB4aorDxDVFaeISWLWIZUWeOIlo4ZUqVj7l810NUqvur00CnAy/cnzRE8/7jYLNuj5ji9LUXk1WJ9u4Gzy5ufyxZO22sBJ5cbDtGH9pBfiVwXdHb8CV7/XsPbCW9yoz/T8VHzPC7j4PL45S8+R8D07YRHrw9vljnAqm74FI+a81O8mON8FP/OV28ucDwZT5Z8pQiKEXQxm6SvKn5VPGEV4285ePRiNp6N57Saz+ejFY1z5LS9XMGbXzf8/8P6xwa+LlZXi9U1jJ7Ct3Xbsg4+EOxqzCKzBUSr7zl8t9j1J63HbjjmJGZncbtp1+3VXXp4AT6Tki84njbNE/40fwH7bqaciw4AAA=='
        _decode_mol_cell('chemical/x-mdl-molfile', 'binary', base64_str)

if __name__ == "__main__":
    main()