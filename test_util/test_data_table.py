"""
Copyright (C) 2017 Anodyne Informatics, LLC
"""

from unittest import TestCase, main

from ruse.util.data_table import DataTable


class TestDataTable(TestCase):
    """A class for testing data table functionality"""

    def test_add_columns(self) -> None:
        """Tests column naming"""

        columns = [DataTable.column_definition('test', 'string'), DataTable.column_definition('test (2)', 'string')]
        table = DataTable(columns=columns, data=[[None, None]])
        table.add_column('test', 'string')
        self.assertEqual(table.columns[2]['name'], 'test (3)')

        columns = [DataTable.column_definition('test (3)', 'string'), DataTable.column_definition('test (2)', 'string')]
        table = DataTable(columns=columns, data=[[None, None]])
        table.add_column('test', 'string')
        self.assertEqual(table.columns[2]['name'], 'test')

if __name__ == "__main__":
    main()

