import concurrent.futures as cf
import itertools
import math
import os
import time

from typing import Any, Optional, Union
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS, rdMolAlign
from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, ColumnData, \
    TableData, string_input_field, integer_input_field, boolean_input_field

from decimal import *

# So that Mol properties, such as the _Name, are pickled when passed
# into the multiprocessing bit.  Needs to be global because of use of
# concurrent.futures.ProcessPoolExecutor.
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


def already_have_mcs(mcss: list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]],
                     new_mcs_dict: dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]) -> bool:
    """
    See if we already have the new_mcs_smarts as a substructure.
    Checking the SMARTS strings isn't enough because they aren't
    canonicalised, so check that all the atom lists match.
    Args:
        mcss ([dict, ]):
        new_mcs_dict (dict):

    Returns:
        bool
    """
    # print(f'Checking if we have {new_mcs_dict["smartsString"]}'
    #       f' {new_mcs_dict["numMols"]} {new_mcs_dict["numAtoms"]}'
    #       f' {new_mcs_dict["numBonds"]}')
    if new_mcs_dict['smartsString'] in mcss:
        return True

    def extract_matches(mcs_dict):
        new_matches = []
        for mol in mcs_dict['mols']:
            matches = mol.GetSubstructMatches(mcs_dict['qmol'])
            match_list = []
            for m in matches:
                match_list.append(sorted(list(m)))
            match_list.sort()
            new_matches.append(match_list)
        new_matches.sort()
        return new_matches

    new_mcs_matches = extract_matches(new_mcs_dict)
    # print(f'new matches : {new_mcs_matches}')

    for mcs in mcss:
        # print(f'Against {mcs["smartsString"]}'
        #       f' {mcs["numMols"]} {mcs["numAtoms"]}'
        #       f' {mcs["numBonds"]}')
        if (new_mcs_dict['numAtoms'] == mcs['numAtoms']
                and new_mcs_dict['numBonds'] == mcs['numBonds']
                and new_mcs_dict['numMols'] == mcs['numMols']):
            these_matches = extract_matches(mcs)
            if these_matches == new_mcs_matches:
                # print(f'XXXXXX {mcs["smartsString"]} and {new_mcs_dict["smartsString"]} the same')
                return True

    return False


def findSingleMCS(sub_mols: list[Chem.Mol], all_mols: list[Chem.Mol],
                  ids: list[Any], min_atom_num: int, min_bond_num: int,
                  complete_rings: bool, start_time: int, max_time: int,
                  threshold: Optional[float] = 1.0) -> tuple[dict[
                                                                 str, Union[str, int, Chem.Mol, list[Chem.Mol]]], bool]:
    """
    Do a single MCS calculation on the given molecules.
    Args:
        sub_mols [Chem.Mol, ]: must be at least 2 molecules in set, obvs.
        all_mols [Chem.Mol, ]: must be at least 2 molecules in set, obvs.
        ids [Any, ]:
        min_atom_num (int): minimum number of atoms in MCS for it to be
                            returned
        min_bond_num (int): minimum number of bonds in MCS for it to be
                            returned
        complete_rings (bool): if True (default) partial ring matches
                              won't be returned.
        start_time (int): when the whole job started
        max_time (int): maximum number of seconds to spend on the thing

    Returns:
        dict of results and bool saying whether it timed out or not.
    """
    no_res = {
        'smartsString': '',
        'numAtoms': 0,
        'numBonds': 0,
        'qmol': None,
    }
    this_max = max_time - int(time.time() - start_time)
    if this_max <= 0:
        return no_res, False
    mcs_res = rdFMCS.FindMCS(sub_mols, completeRingsOnly=complete_rings,
                             threshold=threshold, timeout=this_max)

    if (not mcs_res.smartsString
            or mcs_res.numAtoms < min_atom_num
            or mcs_res.numBonds < min_bond_num):
        return no_res, False
    else:
        mcs_dict = {
            'smartsString': mcs_res.smartsString,
            'numAtoms': mcs_res.numAtoms,
            'numBonds': mcs_res.numBonds,
            'qmol': Chem.MolFromSmarts(mcs_res.smartsString)
        }
        # find any other molecules with this MCS
        hit_mols = []
        hit_ids = []
        if mcs_dict['qmol'] is not None:
            for mol, id in zip(all_mols, ids):
                if mol.HasSubstructMatch(mcs_dict['qmol']):
                    hit_mols.append(mol)
                    hit_ids.append(id)
        mcs_dict['mols'] = hit_mols
        mcs_dict['numMols'] = len(hit_mols)
        mcs_dict['ids'] = hit_ids

    return mcs_dict, mcs_res.canceled


def sort_mcss(mcss: list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]]) -> None:
    """
    Do an in place sort of the MCSS list.
    Order is number of mols that have MCS, number of atoms in MCS,
    number of bonds in MCS.
    Args:
        mcss (list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]]):

    Returns:

    """
    mcss.sort(key=lambda p: p['numBonds'], reverse=True)
    mcss.sort(key=lambda p: p['numAtoms'], reverse=True)
    mcss.sort(key=lambda p: p['numMols'], reverse=True)


def findMCSsExhaustive(mols: list[Chem.Mol], ids: list[Any],
                       min_mol_num: int, min_atom_num: int,
                       min_bond_num: int, complete_rings: Optional[bool] = True,
                       max_time: Optional[int] = 600,
                       num_procs: Optional[int] = -1) -> tuple[
    list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]], bool]:
    """

    Args:
        mols [Chem.Mol, ]: must be at least 2 molecules in set, obvs.
        ids [Any, ]:
        min_mol_num (int): minimum number of molecules having an MCS for
                         it to be returned.
        min_atom_num (int): minimum number of atoms in MCS for it to be
                            returned
        min_bond_num (int): minimum number of bonds in MCS for it to be
                            returned
        complete_rings (bool): if True (default) partial ring matches
                              won't be returned.
        max_time (int): maximum number of seconds to spend on the
                        thing. Defaults to 600s.
        num_procs (int): number of parallel processes to use.  If -ve,
                         uses os.cpu_count() - num_procs.  Default -1.
                         Obviously final count must be at least 1.
    Returns:
        tuple of
        list of {str: dict}: dict of dicts with results of different
                             MCSs found.  Inner dict gives details of
                             each MCS and molecules in it.
        bool: whether it timed out or not.
    """
    mcss = []
    nmols = len(mols)

    def stirling_fac(n):
        twopi = Decimal(2 * math.pi)
        e = Decimal(math.e)
        return Decimal.sqrt(twopi * n) * ((n / e) ** n)

    nmol_fac = stirling_fac(Decimal(nmols))
    min_mol_num_fac = stirling_fac(Decimal(min_mol_num))
    diff_fac = stirling_fac(Decimal(nmols - min_mol_num))
    approx_ncomb = int(nmol_fac / (min_mol_num_fac * diff_fac))
    print(f'Approximate number of MCS determinations : {approx_ncomb}')

    if num_procs < 0:
        num_procs = os.cpu_count() + num_procs
    if num_procs < 1:
        num_procs = 1

    start_time = time.time()
    timed_out = False
    with cf.ProcessPoolExecutor(num_procs) as executor:

        combinator = itertools.combinations(mols, min_mol_num)
        step_size = num_procs * 5
        num_subbed = 0
        while True:
            futures = []
            for i in range(step_size):
                try:
                    comb = next(combinator)
                except StopIteration:
                    break
                if time.time() - start_time > max_time:
                    executor.shutdown(wait=False, cancel_futures=True)
                    timed_out = True
                    break
                sub_mols = list(comb)
                fut = executor.submit(findSingleMCS, sub_mols, mols, ids,
                                      min_atom_num, min_bond_num,
                                      complete_rings, start_time, max_time)
                futures.append(fut)

            num_subbed += len(futures)
            if futures:
                try:
                    for fut in cf.as_completed(futures, max_time):
                        mcs_dict, tim = fut.result()
                        if tim:
                            timed_out = True
                        if mcs_dict['smartsString'] and not already_have_mcs(mcss, mcs_dict):
                            mcss.append(mcs_dict)
                except cf.TimeoutError as err:
                    timed_out = True
                    break
            else:
                break

    sort_mcss(mcss)
    return mcss, timed_out


def findMCSsVGreedy(mols: list[Chem.Mol], ids: list[Any], min_mol_num: int,
                    min_atom_num: int, min_bond_num: int,
                    complete_rings: Optional[bool] = True,
                    max_time: Optional[int] = 600) -> tuple[
    list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]], bool]:
    """

    Args:
        mols [Chem.Mol, ]: must be at least 2 molecules in set, obvs.
        ids [Any, ]: ID from the spotfire column
        min_mol_num (int): minimum number of molecules having an MCS for
                         it to be returned.
        min_atom_num (int): minimum number of atoms in MCS for it to be
                            returned
        min_bond_num (int): minimum number of bonds in MCS for it to be
                            returned
        complete_rings (bool): if True (default) partial ring matches
                              won't be returned.
        max_time (int): maximum number of seconds to spend on the
                        thing. Defaults to 600s.
        num_procs (int): number of parallel processes to use.  If -ve,
                         uses os.cpu_count() - num_procs.  Default -1.
                         Obviously final count must be at least 1.
    Returns:
        tuple of
        list of {str: dict}: dict of dicts with results of different
                             MCSs found.  Inner dict gives details of
                             each MCS and molecules in it.
        bool: whether it timed out or not.
    """
    mcss = []
    start_time = time.time()
    curr_mols = mols
    timed_out = False
    while True:
        thresh = min_mol_num / len(curr_mols)
        mcs, tim = findSingleMCS(curr_mols, mols, ids, min_atom_num, min_bond_num,
                                 complete_rings, start_time, max_time,
                                 threshold=thresh)
        if tim:
            timed_out = True
        if not mcs['smartsString']:
            break
        if not already_have_mcs(mcss, mcs):
            mcss.append(mcs)
        new_curr_mols = [m for m in curr_mols if m not in mcs['mols']]
        if len(new_curr_mols) < min_mol_num:
            break
        curr_mols = new_curr_mols

    sort_mcss(mcss)
    return mcss, timed_out


def findMCSs(mols: list[Chem.Mol], min_mol_num: int, min_atom_num: int,
             min_bond_num: int, ids: list[Any],
             complete_rings: Optional[bool] = True,
             max_time: Optional[int] = 600, num_procs: Optional[int] = -1,
             method: Optional[str] = 'EXHAUSTIVE') -> tuple[
    list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]], bool]:
    """

    Args:
        mols [Chem.Mol, ]: must be at least 2 molecules in set, obvs.
        min_mol_num (int): minimum number of molecules having an MCS for
                         it to be returned.
        min_atom_num (int): minimum number of atoms in MCS for it to be
                            returned
        min_bond_num (int): minimum number of bonds in MCS for it to be
                            returned
        ids ([Any,]):
        complete_rings (bool): if True (default) partial ring matches
                              won't be returned.
        max_time (int): maximum number of seconds to spend on the
                        thing. Defaults to 600s.
        num_procs (int): number of parallel processes to use.  If -ve,
                         uses os.cpu_count() - num_procs.  Default -1.
                         Obviously final count must be at least 1.
        method (str): either EXHAUSTIVE or GREEDY.  The latter is
                      likely to be faster be produce a less complete
                      set of possibilities for the MCSs
    Returns:
        tuple of
        list of {str: dict}: dict of dicts with results of different
                             MCSs found.  Inner dict gives details of
                             each MCS and molecules in it.
        bool: whether it timed out or not.
    """
    if method == 'GREEDY':
        return findMCSsVGreedy(mols, ids, min_mol_num, min_atom_num, min_bond_num,
                               complete_rings, max_time)
    elif method == 'EXHAUSTIVE':
        return findMCSsExhaustive(mols, ids, min_mol_num, min_atom_num, min_bond_num,
                                  complete_rings, max_time, num_procs)
    else:
        print(f'ERROR : method must be either "GREEDY" or "EXHAUSTIVE", you'
              ' gave {method}.')
        return {}, False


def extract_input_values(request: DataFunctionRequest) -> tuple[dict[str, Union[str, int, bool, list[Chem.Mol]]],
                                                                str, DataType]:
    inputs = {}
    column_id = string_input_field(request, 'structureColumn')
    input_column = request.inputColumns[column_id]
    input_column_name = input_column.name
    mols = [m for m in column_to_molecules(input_column) if m is not None]
    inputs['mols'] = mols

    column_id = string_input_field(request, 'idColumn')
    id_column = request.inputColumns[column_id]
    inputs['ids'] = id_column.values
    ids_type = id_column.dataType

    for id in [('minMCSMols', 'min_mol_num'), ('minMCSAtoms', 'min_atom_num'),
               ('minMCSBonds', 'min_bond_num'), ('maxTime', 'max_time')]:
        inputs[id[1]] = integer_input_field(request, id[0])
    inputs['complete_rings'] = not boolean_input_field(request, 'notCompleteRings')
    inputs['method'] = string_input_field(request, 'method')
    return inputs, input_column_name, ids_type


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

    high_bnds = []
    matches = mol.GetSubstructMatches(qmol)
    for match in matches:
        for pair in itertools.combinations(match, 2):
            bond = mol.GetBondBetweenAtoms(pair[0], pair[1])
            if bond is not None:
                high_bnds.append(bond.GetIdx())

    # print(f'atoms : {high_ats}')
    # print(f'bonds : {high_bnds}')
    high_bnds_str = ' '.join([str(b + 1) for b in high_bnds])
    prop_text = f'COLOR #ff0000\nATOMS\nBONDS {high_bnds_str}'
    mol.SetProp('Renderer_Highlight', prop_text)


def build_output_table(input_column_name: str, mcss:
list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]], timed_out: bool,
                       method: str, ids_type: DataType) -> TableData:
    """
    Put together a table from the mcss results.  The rows are arranged
    to show the different MCSs and the molecules that have them.
    Each row has the query mol that defines the MCS, the SMARTS
    equivalent, the number of molecules with that MCS and a molecule
    with that MCS. All members of the MCS are contiguous in the table.
    The MCS is highlighted in each molecule, and the same molecule can
    (for exhaustive search, not for greedy) appear in more than one
    row.
    Args:
        input_column_name (str): 
        mcss (list[dict[str, Union[str, int, Chem.Mol, list[Chem.Mol]]]]):
        timed_out (bool) 

    Returns:
        TableData - the MCS output in table form.
    """
    rdDepictor.SetPreferCoordGen(True)
    if timed_out:
        table_name = f'{input_column_name} PARTIAL MCSs {method}'
    else:
        table_name = f'{input_column_name} MCSs {method}'

    smarts = []
    qmols = []
    num_mols = []
    mols = []
    ids = []
    mcs_nums = []
    mcs_natoms = []
    mcs_nbonds = []
    for i, mcs in enumerate(mcss, 1):
        rdDepictor.Compute2DCoords(mcs['qmol'])
        for mol, id in zip(mcs['mols'], mcs['ids']):
            smarts.append(mcs['smartsString'])
            qmols.append(mcs['qmol'])
            num_mols.append(mcs['numMols'])
            mcs_nums.append(i)
            mcs_natoms.append(mcs['qmol'].GetNumAtoms())
            mcs_nbonds.append(mcs['qmol'].GetNumBonds())
            # take a copy of the molecule, line it up on the MCS qmol
            # and colour the atoms in the MCS
            mol_cp = Chem.Mol(mol)
            rdDepictor.Compute2DCoords(mol_cp)
            match = mol_cp.GetSubstructMatch(mcs['qmol'])
            atom_map = [(m, i) for i, m in enumerate(match)]
            rdMolAlign.AlignMol(mol_cp, mcs['qmol'], atomMap=atom_map)
            highlight_molecule(mol_cp, mcs['qmol'])
            mols.append(mol_cp)
            ids.append(id)

    mcs_num_column = ColumnData(name='MCS Number',
                                dataType=DataType.INTEGER, values=mcs_nums)
    qmol_column = molecules_to_column(qmols, 'MCS', DataType.BINARY)
    smarts_column = ColumnData(name='MCS SMARTS', contentType='chemical/x-smarts',
                               dataType=DataType.STRING, values=smarts)
    num_mols_column = ColumnData(name='MCS Mol. Count', dataType=DataType.INTEGER, values=num_mols)
    natoms_column = ColumnData(name='MCS Num Atoms', dataType=DataType.INTEGER, values=mcs_natoms)
    nbonds_column = ColumnData(name='MCS Num Bonds', dataType=DataType.INTEGER, values=mcs_nbonds)
    mol_column = molecules_to_column(mols, 'MCS Molecule', DataType.BINARY)
    id_column = ColumnData(name='ID MCS Molecule', dataType=ids_type, values=ids)
    all_cols = [mcs_num_column, qmol_column, natoms_column, nbonds_column,
                num_mols_column, smarts_column, mol_column, id_column]
    output_table = TableData(tableName=table_name, columns=all_cols)
    return output_table


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    inputs, input_column_name, ids_type = extract_input_values(request)
    mcss, timed_out = findMCSs(**inputs)

    output_table = build_output_table(input_column_name, mcss, timed_out,
                                      inputs["method"], ids_type)
    response = DataFunctionResponse(outputTables=[output_table])
    return response
