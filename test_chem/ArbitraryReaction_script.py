
from collections import deque
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction

from df.chem_helper import column_to_molecules, \
    molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, \
    string_input_field
from ruse.rdkit.rdkit_utils import string_to_reaction


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
                    Chem.SanitizeMol(new_prod)
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
    rxn_smarts_field = request.inputFields['reactionSmarts']
    rxn_smarts = str(rxn_smarts_field.data)
    if rxn_smarts is not None and rxn_smarts and rxn_smarts != 'None':
        rxn = string_to_reaction(rxn_smarts_field.contentType, rxn_smarts)
    else:
        rxn_sketcher_field = request.inputFields['reactionQuery']
        rxn_sketch_text = str(rxn_sketcher_field.data)
        if rxn_sketch_text is not None and rxn_sketch_text and rxn_sketch_text != 'None':
            rxn = string_to_reaction(rxn_sketcher_field.contentType, rxn_sketch_text)
        else:
            raise ValueError('Received neither reaction SMARTS nor'
                             ' reaction sketch.')

    input_column = request.inputColumns[column_id]
    mols = column_to_molecules(input_column)
    products = run_reactions(mols, rxn)
    products_column = molecules_to_column(products, f'{input_column.name} Products',
                                          DataType.BINARY)
    response = DataFunctionResponse(outputColumns=[products_column])
    return response
