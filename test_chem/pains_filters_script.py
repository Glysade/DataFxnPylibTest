import csv
import itertools
import os
from pathlib import Path
from typing import Optional, Union

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, ColumnData, \
    string_input_field
from rdkit import Chem


def build_query_dict(query_defs: dict[str, str]) -> dict[str, Union[Chem.Mol, str]]:
    """
    Makes RDKit query mols for the SMARTS in query_defs
    Args:
        query_defs [dict]:

    Returns:
        dict
    """
    query_dict = {}
    for name, smt in query_defs.items():
        qmol = Chem.MolFromSmarts(smt)
        if not qmol:
            print(f'ERROR : failed to process SMARTS string {smt} with name {name}.')
        else:
            query_dict[name] = {'qmol': qmol, 'SMARTS': smt}

    return query_dict


def match_queries(mols: list[Chem.Mol], query_dict) -> list[list[str]]:
    """
    Match the queries against the molecules.  Returns list of list of
    strings, 1 entry in outer list for each molecule, giving the
    names of queries that matched.
    Args:
        mols ([Chem.Mol, ]):
        query_dict ({}):

    Returns:
        [[str,], ]
    """
    hits = []
    for i, mol in enumerate(mols):
        mol_hits = []
        if mol is not None and mol:
            for name, query in query_dict.items():
                if mol.HasSubstructMatch(query['qmol']):
                    mol_hits.append(name)
                    # try:
                    #     print(f'{mol.GetProp("_Name")} : {name}')
                    # except KeyError:
                    #     print(f'Molecule number {i} : {name}')
        hits.append(mol_hits)

    return hits


def read_pains_queries() -> Optional[dict[str, str]]:
    """
    Read the PAINS SMARTS patterns from
    $RDBASE/Data/Pains/wehi_pains.csv.  Returns None if that file not
    found.
    Returns:
        Dict[str: str]: the SMARTS keyed on name.
    """
    try:
        rdbase = os.environ['RDBASE']
    except KeyError:
        print(f'ERROR : no RDBASE')
        return None

    print(os.environ['RDBASE'])
    pains_defs_file = Path(rdbase) / 'Data' / 'Pains' / 'wehi_pains.csv'
    queries = {}
    try:
        with open(pains_defs_file, 'r', newline='') as f:
            csvreader = csv.reader(f)
            for row in csvreader:
                smt_name = row[1].replace('<regId=', '').replace('>', '')
                queries[smt_name] = row[0]

    except IOError:
        print(f'ERROR : no file {pains_defs_file}')
        return None
    return queries


def highlight_molecule(mol: Chem.Mol, qmol: Chem.Mol) -> None:
    """
    Use the qmol to add a highlight tag to mol in the Glysade format
    for display in Spotfire.
    Args:
        mol (Chem.Mol): the molecule to be highlighted.
        smarts (str): the query that defines the atoms and bonds to
                         be highlighted.

    Returns:

    """
    if not qmol or qmol is None:
        return

    high_ats = []
    high_bnds = []
    # If multiple PAINS hit, and the hits overlap on a molecule, they
    # don't show up unless the individual SMARTS patterns are matched
    # separately.
    for qmol_frag in Chem.GetMolFrags(qmol, asMols=True):
        matches = mol.GetSubstructMatches(qmol_frag)
        for match in matches:
            for pair in itertools.combinations(match, 2):
                bond = mol.GetBondBetweenAtoms(pair[0], pair[1])
                if bond is not None:
                    high_bnds.append(bond.GetIdx())
            high_ats.extend(match)

    high_ats = list(set(high_ats))
    high_bnds = list(set(high_bnds))
    # print(f'atoms : {high_ats}')
    # print(f'bonds : {high_bnds}')
    high_ats_str = ' '.join([str(a + 1) for a in high_ats])
    high_bnds_str = ' '.join([str(b + 1) for b in high_bnds])
    prop_text = f'COLOR #ff0000\nATOMS {high_ats_str}\nBONDS {high_bnds_str}'
    mol.SetProp('Renderer_Highlight', prop_text)


def highlight_molecules(mols: list[Chem.Mol], all_smarts: list[Chem.Mol]) -> None:
    for mol, smarts in zip(mols, all_smarts):
        if mol is not None:
            highlight_molecule(mol, smarts)


def run_pains(mols: list[Chem.Mol]) -> list[tuple[bool, str, str]]:
    """
    Run the PAINS filters on each molecule.  Returns a list of tuples.
    Each tuple has a bool for whether it had a PAINS result,
    a string which is the comma-separated list of names of any hits and
    a string which is the dot-separated list of SMARTS patterns of any
    hits.
    Args:
        mols ([Chem.Mol, ]):

    Returns:
        [(bool, str, str), ]
    """
    pains_smarts = read_pains_queries()
    pains_queries = build_query_dict(pains_smarts)
    p_hits = match_queries(mols, pains_queries)
    results = []
    for mol, p_h in zip(mols, p_hits):
        res_str = []
        res_smt = []
        res_bool = False
        if p_h:
            res_str.extend(p_h)
            res_smt.extend([pains_queries[p]["SMARTS"] for p in p_h])
            res_bool = True
        all_smt = '.'.join(res_smt)
        if all_smt:
            pains_mol = Chem.MolFromSmarts(all_smt)
        else:
            pains_mol = None
        results.append((res_bool, ','.join(res_str), pains_mol))

    return results


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'structureColumn')
    input_column = request.inputColumns[column_id]
    mols = column_to_molecules(input_column)
    pains = run_pains(mols)
    pains_cols = [p[0] for p in pains]
    pains_column = ColumnData(name=f'{input_column.name} PAINS?', dataType=DataType.BOOLEAN,
                              values=pains_cols)
    pains_strs = [p[1] for p in pains]
    pains_strs_column = ColumnData(name=f'{input_column.name} PAINS NAMES', dataType=DataType.STRING,
                                   values=pains_strs)
    pains_mols = [p[2] for p in pains]
    pains_smarts_column = molecules_to_column(pains_mols, f'{input_column.name} PAINS SMARTS',
                                              DataType.BINARY)
    highlight_molecules(mols, pains_mols)
    pains_hits_column = molecules_to_column(mols, f'{input_column.name} PAINS MAP', DataType.BINARY)
    response = DataFunctionResponse(outputColumns=[pains_column, pains_strs_column,
                                                   pains_smarts_column, pains_hits_column])
    return response
