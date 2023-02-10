import os
from unittest import TestCase, main

from rdkit import Chem
from rdkit.Chem import AllChem

from ruse.rdkit.mmp import transform, result_to_mmp_transform, environment_rule_pairs, create_combined_structure, \
    align_combined_molecules, database_property_names
from ruse.rdkit.rdkit_utils import remove_atom_mappings, remove_explicit_hydrogens, RDKitFormat

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def default_database_filepath() -> str:
    return database_filepath('bioactivity-18_22_53_07')


def large_database_filepath() -> str:
    return database_filepath('AID_2642_datatable_all')


def database_with_missing() -> str:
    return database_filepath('AID_1341_datatable_all')


def database_filepath(base: str) -> str:
    if os.name == 'nt':
        return os.path.join('c:/', 'db', 'mmp', f'{base}.mmpdb')
    return os.path.join(os.sep, 'mnt', 'c', 'db', 'mmp', f'{base}.mmpdb')


class TestMmp(TestCase):

    def test_transform(self) -> None:
        database_file = default_database_filepath()
        smiles = 'c1cccnc1O'
        result = transform(smiles, ['ACTIVITY'], database_file)
        cnt = 0
        for product in result.transform_products:
            for rule in product.property_rules:
                if not rule:
                    continue
                cnt += 1
                self.eval_result(smiles, product, rule)
        self.assertEqual(cnt, 263)

    def test_transform2(self) -> None:
        database_file = default_database_filepath()
        smiles = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)Cc5ccc(O)cc5)nc43)C(O)C2O)n1'
        result = transform(smiles, ['ACTIVITY'], database_file)
        cnt = 0
        for product in result.transform_products:
            for rule in product.property_rules:
                if not rule:
                    continue
                cnt += 1
                if False:
                    # information for debugging
                    print('{} {} {} {} >> {} {} {}| {} {} {}'.format(cnt, product.smiles, rule.radius, rule.from_smiles,
                                                                     rule.to_smiles,
                                                                     rule.avg, rule.is_reversed, rule.constant_smiles,
                                                                     rule.variable_smiles,
                                                                     rule.attachment_order))
                    if len(rule.attachment_order) > 2 or cnt == 421:
                        print('structure = "{}"'.format(smiles))
                        print('product = "{}"'.format(product.smiles))
                        print('from_smiles = "{}"'.format(rule.from_smiles))
                        print('to_smiles = "{}"'.format(rule.to_smiles))
                        print('reversed = True if 1 == {} else False'.format(rule.is_reversed))
                        print('constant_smiles = "{}"'.format(rule.constant_smiles))
                        print('variable_smiles = "{}"'.format(rule.variable_smiles))
                        print('attachment_order = [int(c) for c in "{}"]'.format(rule.attachment_order))
                self.eval_result(smiles, product, rule)

        self.assertEqual(cnt, 1274)

    def assert_equal_smiles(self, reference_mol, combination_mol):
        reference_smiles = Chem.MolToSmiles(reference_mol, True)
        combination_smiles = Chem.MolToSmiles(remove_atom_mappings(combination_mol), True)
        test = reference_smiles == combination_smiles
        # because of the symmetry of the 2 stereo centers these guys are the same, but I haven't figured out how to
        # code this up in RDKit (see the notebook mmp_create_structure.ipynb)
        if not test:
            if reference_smiles == 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)C[C@]56CC[C@](O)(CC5)CC6)nc43)C(O)C2O)n1' \
                    and combination_smiles == 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)C[C@@]56CC[C@@](O)(CC5)CC6)nc43)C(O)C2O)n1':
                return
            reference_smiles = Chem.MolToSmiles(remove_explicit_hydrogens(reference_mol), True)
            test = reference_smiles == combination_smiles
        self.assertTrue(test)

    def eval_result(self, structure, product, rule):
        product = product.smiles
        to_smiles = rule.to_smiles
        constant_smiles = rule.constant_smiles
        variable_smiles = rule.variable_smiles
        combined_mol = create_combined_structure(constant_smiles, variable_smiles)
        combined_prod = create_combined_structure(constant_smiles, to_smiles)
        structure_mol = Chem.MolFromSmiles(structure)
        product_mol = Chem.MolFromSmiles(product)
        AllChem.Compute2DCoords(structure_mol)
        align_combined_molecules(structure_mol, combined_mol, product_mol, constant_smiles)

        self.assert_equal_smiles(structure_mol, combined_mol)
        self.assert_equal_smiles(product_mol, combined_prod)

    def test_trans_bond_in_variable(self):
        # this test in not applicable to MMPDB2 as we don't need to include the attachment order.
        structure = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)Cc5ccc(O)cc5)nc43)C(O)C2O)n1'
        product = 'CCN/C(=N\CO)NCCCN(Cc1ccc(O)cc1)c1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1'
        to_smiles = '[*:1]/N=C(\\NCC)NCCCN([*:2])[*:3]'
        constant_smiles = '[*]CO.[*]Cc1ccc(O)cc1.[*]c1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1'
        variable_smiles = '[*]NC([*])[*]'
        attachment_order = [int(c) for c in "201"]

        combined_prod = create_combined_structure(constant_smiles, to_smiles)
        combined_mol = create_combined_structure(constant_smiles, variable_smiles, attachment_order)
        structure_mol = Chem.MolFromSmiles(structure)
        product_mol = Chem.MolFromSmiles(product)

        structure_smiles = Chem.MolToSmiles(remove_atom_mappings(structure_mol), True)
        combined_smiles = Chem.MolToSmiles(remove_atom_mappings(combined_mol), True)
        product_smiles = Chem.MolToSmiles(remove_atom_mappings(product_mol), True)
        combined_prod_smiles = Chem.MolToSmiles(remove_atom_mappings(combined_prod), True)

        self.assertEqual(structure_smiles, combined_smiles)
        self.assertEqual(product_smiles, combined_prod_smiles)

    def test_to_smiles_ordering(self):
        product = "CCc1cc(O)ccc1C(CO)Nc1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1"
        to_smiles = "[*:2]c1ccc([*:1])c(CC)c1"
        constant_smiles = "[*]C(CO)Nc1nc(N)c2ncn(C3OC(c4nnn(CC)n4)C(O)C3O)c2n1.[*]O"

        combined_prod = create_combined_structure(constant_smiles, to_smiles)
        product_mol = Chem.MolFromSmiles(product)
        combined_prod = remove_atom_mappings(combined_prod)

        product_smiles = Chem.MolToSmiles(product_mol, True)
        combined_prod_smiles = Chem.MolToSmiles(combined_prod, True)
        self.assertEqual(product_smiles, combined_prod_smiles)

    def test_trans_bond_in_variable2(self):
        structure = "CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)Cc5ccc(O)cc5)nc43)C(O)C2O)n1"
        product = "CCn1nnc(C2OC(n3cnc4c(/C(N)=N/O)nc(NC(CO)Cc5ccc(O)cc5)nc43)C(O)C2O)n1"
        from_smiles = "[*:1]N"
        to_smiles = "[*:1]/C(N)=N/O"
        reversed = True if 1 == 1 else False
        constant_smiles = "[*]c1nc(NC(CO)Cc2ccc(O)cc2)nc2c1ncn2C1OC(c2nnn(CC)n2)C(O)C1O"
        variable_smiles = "[*]N"

        combined_prod = create_combined_structure(constant_smiles, to_smiles)
        combined_mol = create_combined_structure(constant_smiles, variable_smiles)
        structure_mol = Chem.MolFromSmiles(structure)
        product_mol = Chem.MolFromSmiles(product)

        structure_smiles = Chem.MolToSmiles(remove_atom_mappings(structure_mol), True)
        combined_smiles = Chem.MolToSmiles(remove_atom_mappings(combined_mol), True)
        product_smiles = Chem.MolToSmiles(remove_atom_mappings(product_mol), True)
        combined_prod_smiles = Chem.MolToSmiles(remove_atom_mappings(combined_prod), True)

        self.assertEqual(structure_smiles, combined_smiles)
        self.assertEqual(product_smiles, combined_prod_smiles)

    def test_create_mmp_transform1(self):
        self.eval_create_mmp_transform('c1cccnc1O', 263)

    def test_create_mmp_transform2(self):
        self.eval_create_mmp_transform('CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)Cc5ccc(O)cc5)nc43)C(O)C2O)n1', 1295)

    def eval_create_mmp_transform(self, query: str, expected_product_count: int) -> None:
        properties = ['PSA', 'ACTIVITY']
        database_file = default_database_filepath()
        result = transform(query, properties, database_file)
        query_mol = Chem.MolFromSmiles(query)
        transforms = result_to_mmp_transform(query_mol, result)
        grid = transforms.to_grid(molecules=True, column_major=True)
        self.assertEqual(8, len(grid))
        self.assertEqual(expected_product_count, len(grid[0]))
        data_table = transforms.to_data_table(RDKitFormat.smi, query)
        self.assertEqual(len(data_table.data), expected_product_count)
        self.assertEqual(10, len(data_table.data[0]))
        self.assertEqual(10, len(data_table.columns))

        for product in transforms.products:
            for property, tp in zip(properties, product.transform_products):
                environment = environment_rule_pairs(tp.rule_environment_id, tp.reversed, property, database_file)
                delta1 = tp.value
                delta2 = environment.delta()
                self.assertAlmostEqual(delta1, delta2)

    def test_transform_large_database(self):
        database_file = large_database_filepath()
        self.assertTrue(os.path.isfile(database_file), 'MMPDB large test database {} not built'.format(database_file))
        properties = database_property_names(database_file)
        self.assertEqual(set(properties), {'ACTIVITY', 'MOLWEIGHT', 'LOGP', 'TPSA'})
        query = 'CC(C)CN(C1=C(N(C(=O)NC1=O)CC2=CC=CC=C2)N)C(=O)CCC3=CC=CC=C3Cl'
        result = transform(query, properties, database_file)
        actual = len(result.transform_products)
        expected = 14450
        self.assertTrue(expected - 1000 <= actual <= expected + 1000)

    def test_large_database1(self) -> None:
        query = 'CC(C)CN(C1=C(N(C(=O)NC1=O)CC2=CC=CC=C2)N)C(=O)CCC3=CC=CC=C3Cl'
        self.eval_large_database(query, 14500, 14500)

    def test_large_database2(self) -> None:
        query = 'CCn1nnc(C2OC(n3cnc4c(N)nc(NC(CO)Cc5ccc(O)cc5)nc43)C(O)C2O)n1'
        self.eval_large_database(query, 9500, 9500)

    def eval_large_database(self, query, expected_product_count, expected_transform_count):
        properties = ['TPSA', 'ACTIVITY']
        database_file = large_database_filepath()
        self.assertTrue(os.path.isfile(database_file), 'MMPDB large test database {} not built'.format(database_file))
        result = transform(query, properties, database_file)
        query_mol = Chem.MolFromSmiles(query)
        transforms = result_to_mmp_transform(query_mol, result)
        grid = transforms.to_grid(molecules=True, column_major=True)
        data_table = transforms.to_data_table(RDKitFormat.smi, query)
        # print('product count {} transform count {} expected product count {} expected transform count {}'.format(
        #     len(grid[0]), len(transforms.products), expected_product_count, expected_transform_count))
        self.assertTrue(expected_transform_count - 1000 <= len(transforms.products) <= expected_transform_count + 1000)
        self.assertEqual(len(grid), 8)
        self.assertTrue(expected_product_count - 1000 <= len(grid[0]) <= expected_product_count + 1000)
        self.assertEqual(len(data_table.data), len(grid[0]))
        self.assertEqual(len(data_table.data[0]), 10)
        self.assertEqual(len(data_table.columns), 10)

        for product in transforms.products:
            for property, tp in zip(properties, product.transform_products):
                environment = environment_rule_pairs(tp.rule_environment_id, tp.reversed, property, database_file)
                delta1 = tp.value
                delta2 = environment.delta()
                self.assertAlmostEqual(delta1, delta2)

    def test_database_with_missing(self):
        database_file = database_with_missing()
        self.assertTrue(os.path.isfile(database_file), 'MMPDB database with missing {} not built'.format(database_file))
        properties = ['ACTIVITY', 'MOLWEIGHT']
        query = 'C1=CC=C2C(=C1)N=C3C(=C(N(C3=N2)/N=C/C4=CC(=C(C=C4)O)O)N)S(=O)(=O)C5=CC(=CC=C5)Cl'
        result = transform(query, properties, database_file)
        query_mol = Chem.MolFromSmiles(query)
        transforms = result_to_mmp_transform(query_mol, result)
        grid = transforms.to_grid(molecules=True, column_major=True)
        data_table = transforms.to_data_table(RDKitFormat.smi, query)
        self.assertEqual(10, len(data_table.data[0]))
        self.assertEqual(114, len(data_table.data))

        properties = ['ACTIVITY']
        result = transform(query, properties, database_file)
        query_mol = Chem.MolFromSmiles(query)
        transforms = result_to_mmp_transform(query_mol, result)
        grid = transforms.to_grid(molecules=True, column_major=True)
        data_table = transforms.to_data_table(RDKitFormat.smi, query)
        self.assertEqual(9, len(data_table.data[0]))
        self.assertEqual(57, len(data_table.data))

    def test_property_names(self) -> None:
        database_file = default_database_filepath()
        property_names = database_property_names(database_file)
        self.assertEqual(set(property_names), {'PSA', 'ACTIVITY', 'MOLWEIGHT', 'ALOGP'})


if __name__ == "__main__":
    main()
