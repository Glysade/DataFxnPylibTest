import re

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, ColumnData, \
    string_input_field
from rdkit import Chem

# These SMARTS are taken from
# https://github.com/OpenEye-Contrib/SmiV/blob/master/test_dir/hbond.smt
# and are subject to the license at
# https://github.com/OpenEye-Contrib/SmiV/blob/master/LICENSE
# viz:
'''
Unless otherwise noted, all files in this directory and all
subdirectories are distributed under the following license:

Copyright (C) 2016
AstraZeneca, David Cosgrove

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met: 

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following 
      disclaimer in the documentation and/or other materials provided 
      with the distribution.
    * Neither the name of AstraZeneca nor the names of its 
      contributors may be used to endorse or promote products derived 
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


The SMILES files in the test_dir directory were extracted from Chembl
(www.ebi.ac.uk/chembl).  None of them are subject to this license.
'''
# The SMARTS are the same as in the original file, the format
# has been altered slightly for ease of use in Python.
HBOND_SMARTS = '''#
#  Definitions of h-bond donors and acceptors
#  Originally by Pete Kenny, probably.
#
#  Hydrogen bond donors and acceptors are defined first
#  without reference to charge.  Charged donors and 
#  acceptors are defined for use in their own right and to
#  specify neutral acceptors and donors. Ylids are defined
#  as a special case because CORINA converts pentavalent
#  representation of hypervalent nitrogen to ylid in the
#  generation of 3D coordinates. 
#
#----------------------------------------------------------
#  Define ylids beacuse CORINA will convert pentavalent
#  nitrogen form of hypervalent species to ylid form
Ylid    [A-]~[N,n;+;H0]
#----------------------------------------------------------
# Define donors (neutral & cationic) 
AllDon  [O,N,n;!H0]
#----------------------------------------------------------
#  Define acceptors (neutral & anionic)
Ac1     [O;X1]
Ac2     [OH]
Ac3     O([CX4])[CX4]
OAc     [$Ac1,$Ac2,$Ac3]
Ac4     [nX2]
Ac5     [NX2]=C
Ac6     [N;X3;H2][CX4]
Ac7     [N;X3;H]([CX4])[CX4]
Ac8     [N;X3]([CX4])([CX4])[CX4]
Ac9     [N;X1]
# Adding in SMARTS to take out none acceptors in oxadiazoles
NonAccN  [$(ncC=O)]
NAc     [$Ac4,$Ac5,$Ac6,$Ac7,$Ac8,$Ac9;!$NonAccN]
AllAcc  [$OAc,$NAc]
#----------------------------------------------------------
# Define catonic donors
Po1     [N,n;+;H,H2,H3]
Po2     [A;H,H2]-C=[$Po1]
PosDon  [$Po1,$Po2]
#----------------------------------------------------------
# Define anionic acceptors
Ne1     [*;-]
Ne2     [O,S;X1]=*[$Ne1]
Ne3     [nX2][n-]
Ne4     [nX2][nX2][n-]
Ne5     [nX2][nX2][nX2][n-]
Ne6     [O,S;X1]=[$(ccc[A-]),$(CC=C[A-])]
NegAcc  [$Ne1,$Ne2,$Ne3,$Ne4,$Ne5,$Ne6;!$Ylid]
#----------------------------------------------------------
# Define neutral donors & acceptors by eliminating 
# charged species from all donors & acceptors
NeuDon  [$AllDon;!$PosDon]
NeuAcc  [$AllAcc;!$NegAcc]
'''
def expand_smarts(smt: str, smarts_dict: dict[str: str]) -> str:
    """
    Take the SMARTS pattern and replace any vector bindings with the
    corresponding SMARTS pattern until there are none left.
    The Daylight theory manual says a vector binding name matches this regex:
    /[a-zA-Z][a-zA-Z0-9_]*/.
    If there's an error, it returns the original SMARTS and possibly
    a second string with an error message or None if it parse ok.
    Args:
        smt (str): SMARTS string to be expanded

    Returns:
        str: expanded SMARTS string
    """

    VBP = re.compile(r'\$[a-zA-Z][a-zA-Z0-9_]*')
    while True:
        m = re.search(VBP, smt)
        if not m:
            break
        vb = smt[m.start() + 1:m.end()]
        smt = smt[:m.start()] + '$(' + smarts_dict[vb] + ')' + smt[m.end():]

    return smt


def read_hbond_queries() -> dict[str: str]:
    """
    Read the Hbonds SMARTS patterns from HBONDS_SMARTS above.
    Returns:
        Dict[str: str]: the SMARTS keyed on name.
    """
    queries = {}
    multi_space = re.compile(r'\s+')
    for smarts_line in HBOND_SMARTS.split('\n'):
        if not smarts_line or smarts_line.startswith('#'):
            continue
        (smt_name, smt_def) = multi_space.split(smarts_line.strip())
        queries[smt_name] = smt_def

    return queries


def highlight_atoms_in_molecule(mol: Chem.Mol, atoms1: list[int],
                                atoms2: list[int]) -> None:
    """
    Highlight just the atoms in the two lists.  Because they are
    donors and acceptors, they will all be individual atoms so
    no bonds need highlighting.  Atoms in atoms1 will be red,
    atoms2 will be blue and atoms in both lists (donors and acceptors
    such as the atom in a hydroxyl), will be orange.
    """
    # The highlighting code counts from 1, whereas the atom lists
    # come in counting from 0.
    red_atoms = set(atoms1)
    blue_atoms = set(atoms2)
    both_atoms = red_atoms.intersection(blue_atoms)
    red_atoms = red_atoms.difference(both_atoms)
    blue_atoms = blue_atoms.difference(both_atoms)

    red_list = ' '.join([str(a+1) for a in red_atoms])
    blue_list = ' '.join([str(a + 1) for a in blue_atoms])
    orange_list = ' '.join([str(a + 1) for a in both_atoms])

    high_str = ''
    if red_list:
        high_str += f'COLOR #dc143c\nATOMS {red_list}\nBONDS\n'
    if blue_list:
        high_str += f'COLOR #00bfff\nATOMS {blue_list}\nBONDS\n'
    if orange_list:
        high_str += f'COLOR #ffbf00\nATOMS {orange_list}\nBONDS\n'
    if high_str:
        mol.SetProp('Renderer_Highlight', high_str)


def find_matches(mol: Chem.Mol, query: Chem.Mol) -> list[int]:
    matches = mol.GetSubstructMatches(query)
    hit_atoms = [n for m in matches for n in m]
    return hit_atoms


def run_hbonds_search(mols: list[Chem.Mol]) -> list[tuple[Chem.Mol, int, int]]:
    """
    Match the hbonds SMARTS against the molecules, returning the counts
    of donor and acceptor matches and a copy of the molecule with the
    Renderer_Highlight prop to show the donors in blue and the
    acceptors in red.
    """
    hbond_queries = read_hbond_queries()
    donq = Chem.MolFromSmarts(hbond_queries['AllDon'])
    accq = Chem.MolFromSmarts(expand_smarts(hbond_queries['AllAcc'],
                                            hbond_queries))
    hbond_results = []
    for mol in mols:
        if mol is not None:
            don_atoms = find_matches(mol, donq)
            acc_atoms = find_matches(mol, accq)
            mol_cp = Chem.Mol(mol)
            highlight_atoms_in_molecule(mol_cp, don_atoms, acc_atoms)
            hbond_results.append((mol_cp, len(don_atoms), len(acc_atoms)))
        else:
            hbond_results.append((None, 0, 0))

    return hbond_results


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    column_id = string_input_field(request, 'structureColumn')
    input_column = request.inputColumns[column_id]
    mols = column_to_molecules(input_column)
    hbonds_counts = run_hbonds_search(mols)
    out_mols = []
    don_counts = []
    acc_counts = []
    for hbc in hbonds_counts:
        out_mols.append(hbc[0])
        don_counts.append(hbc[1])
        acc_counts.append(hbc[2])
    out_mols_col = molecules_to_column(out_mols,
                                       f'Hbonds {input_column.name}',
                                       DataType.BINARY)
    don_counts_col = ColumnData(name=f'Donor counts {input_column.name}',
                                dataType=DataType.INTEGER,
                                values=don_counts)
    acc_counts_col = ColumnData(name=f'Acceptor counts {input_column.name}',
                                dataType=DataType.INTEGER,
                                values=acc_counts)
    response = DataFunctionResponse(outputColumns=[out_mols_col,
                                                   don_counts_col,
                                                   acc_counts_col])

    return response
