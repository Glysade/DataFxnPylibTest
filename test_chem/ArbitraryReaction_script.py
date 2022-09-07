
from collections import deque
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import KekulizeException, AtomValenceException, MolSanitizeException

from df.chem_helper import column_to_molecules, \
    molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, \
    string_input_field
from ruse.rdkit.rdkit_utils import sanitize_mol, string_to_reaction


def highlight_product(mol: Chem.Mol) -> None:
    """
    Takes the molecule, and highlights any atoms that have an
    'old_mapno' property.
    """
    high_ats = []
    for atom in mol.GetAtoms():
        if atom.HasProp('old_mapno'):
            high_ats.append(atom.GetIdx())

    high_bnds = []
    for at1 in high_ats:
        for at2 in high_ats:
            if at1 > at2:
                bond = mol.GetBondBetweenAtoms(at1, at2)
                if bond is not None:
                    high_bnds.append(bond.GetIdx())
    high_ats_str = ' '.join([str(a + 1) for a in high_ats])
    high_bnds_str = ' '.join([str(b + 1) for b in high_bnds])
    prop_text = f'COLOR #ff0000\nATOMS {high_ats_str}\nBONDS {high_bnds_str}'
    mol.SetProp('Renderer_Highlight', prop_text)


def run_reactions(mols: list[Chem.Mol], rxn: ChemicalReaction) -> list[Optional[Chem.Mol]]:
    """
    Run the reaction on each of the input molecules, returning the
    products. The product molecule may be more than 1 fragment, if the
    reaction cleaved a bond, for example.  The reaction may match more
    than one place in a molecule.  In these cases, all reactions are
    applied to the same molecule, so that there will only be 1 product
    for each input molecule.  If there is no product, None is returned
    in the list.

    :param mols ([Chem.Mol, ]:
    :param rxn [ChemicalReaction]:
    :return [[Chem.Mol, ], ]:
    """
    all_prods = []
    for mol in mols:
        if mol is not None and mol:
            final_mols = deque([mol])
            final_mols_smiles = deque([Chem.MolToSmiles(mol)])
            while True:
                next_mol = final_mols.popleft()
                these_prods = rxn.RunReactants((next_mol,))
                # these_prods is a tuple of tuples of Chem.Mol
                for prods in these_prods:
                    new_prod = Chem.Mol()
                    for prod in prods:
                        new_prod = Chem.CombineMols(new_prod, prod)
                    try:
                        Chem.SanitizeMol(new_prod)
                    except MolSanitizeException:
                        continue
                    # sometimes we'll get a molecule that sanitizes but
                    # generates a bad SMILES.  Skip those.
                    smi = Chem.MolToSmiles(new_prod)
                    test_mol = Chem.MolFromSmiles(smi)
                    if test_mol is None:
                        continue
                    highlight_product(new_prod)
                    prod_smi = Chem.MolToSmiles(new_prod)
                    if prod_smi not in final_mols_smiles:
                        final_mols.append(new_prod)
                        final_mols_smiles.append(prod_smi)
                if len(final_mols) < 2:
                    break
            if final_mols:
                all_prods.append(final_mols[0])
            else:
                all_prods.append(None)
        else:
            all_prods.append(None)

    return all_prods


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'structureColumn')
    rxn_sketcher_field = request.inputFields['reactionQuery']
    rxn_sketch_text = str(rxn_sketcher_field.data)
    rxn = string_to_reaction(rxn_sketcher_field.contentType, rxn_sketch_text)

    input_column = request.inputColumns[column_id]
    mols = column_to_molecules(input_column)
    products = run_reactions(mols, rxn)
    products_column = molecules_to_column(products, f'{input_column.name} Products',
                                          DataType.BINARY)
    response = DataFunctionResponse(outputColumns=[products_column])
    return response
