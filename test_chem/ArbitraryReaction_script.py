
from collections import defaultdict, deque
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import MolSanitizeException

from df.chem_helper import column_to_molecules, \
    molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, \
    TableData, string_input_field
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


def combine_prods_to_one_mol(prods: tuple[Chem.Mol]) -> Optional[Chem.Mol]:
    new_prod = Chem.Mol()
    for prod in prods:
        new_prod = Chem.CombineMols(new_prod, prod)
    try:
        Chem.SanitizeMol(new_prod)
    except MolSanitizeException:
        return None
    # sometimes we'll get a molecule that sanitizes but
    # generates a bad SMILES.  Skip those.
    smi = Chem.MolToSmiles(new_prod)
    return Chem.MolFromSmiles(smi)


def run_reaction_single_product(mol: Chem.Mol, rxn: ChemicalReaction,
                                reaction_mode: str) -> list[Optional[Chem.Mol]]:
    """
    Do singleReaction or exhaustiveReaction mode
    """
    final_mols = deque([mol])
    final_mols_smiles = deque([Chem.MolToSmiles(mol)])
    while True:
        next_mol = final_mols.popleft()
        these_prods = rxn.RunReactants((next_mol,))
        # these_prods is a tuple of tuples of Chem.Mol
        for prods in these_prods:
            new_prod = combine_prods_to_one_mol(prods)
            if new_prod is None:
                continue
            prod_smi = Chem.MolToSmiles(new_prod)
            if prod_smi not in final_mols_smiles:
                highlight_product(new_prod)
                final_mols.append(new_prod)
                final_mols_smiles.append(prod_smi)
                if reaction_mode == 'singleReaction':
                    break
        if len(final_mols) < 2:
            break
    if final_mols:
        return [final_mols[0]]
    else:
        return [None]


def run_reaction_multiple_products(mol: Chem.Mol, rxn: ChemicalReaction) -> list[Optional[Chem.Mol]]:
    """
    Do multipleReaction mode
    """
    all_prods = []
    prod_smis = defaultdict(int)
    next_bases = [Chem.Mol(mol)]
    while True:
        next_prods = []
        for next_base in next_bases:
            these_prods = rxn.RunReactants((next_base,))
            for prods in these_prods:
                new_prod = combine_prods_to_one_mol(prods)
                if new_prod is None:
                    continue
                new_prod_smi = Chem.MolToSmiles(new_prod)
                if new_prod_smi not in prod_smis:
                    all_prods.append(new_prod)
                    prod_smis[new_prod_smi] = 1
                    next_prods.append(new_prod)
        if not next_prods:
            break
        next_bases = next_prods

    return all_prods


def run_reaction_on_mol(mol: Chem.Mol, rxn: ChemicalReaction,
                        reaction_mode: str) -> list[Optional[Chem.Mol]]:
    """
    Perform the reaction on a single molecule in the given mode.
    See run_reactions for the explanation of the modes
    """
    if reaction_mode == 'multipleReaction':
        return run_reaction_multiple_products(mol, rxn)
    else:
        return run_reaction_single_product(mol, rxn, reaction_mode)


def run_reactions(mols: list[Chem.Mol], rxn: ChemicalReaction,
                  reaction_mode: str) -> list[tuple[Chem.Mol, Optional[Chem.Mol]]]:
    """
    Run the reaction on each of the input molecules, returning the
    products. The product molecule may be more than 1 fragment, if the
    reaction cleaved a bond, for example. If there is no product, None
    is returned in the list.
    If there is more than 1 reaction site in the molecule, the products
    depends on reaction_mode.
    If reaction_mode is 'singleReaction' then just the first product is
    returned.
    If reaction_mode is 'exhaustiveReaction' then all reactions are
    performed on the single molecule, to return only 1 product.
    If reaction_mode is 'multipleReaction' all reactions are performed
    on successive products to produce multiple products.

    :param mols ([Chem.Mol, ]:
    :param rxn [ChemicalReaction]:
    :param reaction_mode [str]:
    :return [[Chem.Mol, ], ]:
    """
    all_prods = []
    for mol in mols:
        if mol is not None and mol:
            prod_mols = run_reaction_on_mol(mol, rxn, reaction_mode)
            if prod_mols:
                pms = [(mol, pm) for pm in prod_mols]
                all_prods.extend(pms)
            else:
                all_prods.append((mol, None))
        else:
            all_prods.append((mol,None))

    return all_prods


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'structureColumn')
    rxn_sketcher_field = request.inputFields['reactionQuery']
    rxn_sketch_text = str(rxn_sketcher_field.data)
    rxn = string_to_reaction(rxn_sketcher_field.contentType, rxn_sketch_text)

    input_column = request.inputColumns[column_id]
    mols = column_to_molecules(input_column)
    reaction_mode = string_input_field(request, 'reactionMode')

    products = run_reactions(mols, rxn, reaction_mode)
    if reaction_mode == 'multipleReaction':
        init_mols = [p[0] for p in products]
        prods = [p[1] for p in products]
        init_col = molecules_to_column(init_mols, f'{input_column.name}', DataType.BINARY)
        prods_col = molecules_to_column(prods, f'{input_column.name} Products', DataType.BINARY)
        output_table = TableData(tableName=f'{input_column.name} Products',
                                 columns=[init_col, prods_col])
        response = DataFunctionResponse(outputTables=[output_table])
    else:
        prods = [p[1] for p in products]
        products_column = molecules_to_column(prods, f'{input_column.name} Products',
                                              DataType.BINARY)
        response = DataFunctionResponse(outputColumns=[products_column])
    return response
