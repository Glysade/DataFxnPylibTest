from typing import Union

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize
from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunctionRequest, DataFunctionResponse, DataType, \
    string_input_field, boolean_input_field
from ruse.rdkit.rdkit_utils import sanitize_mol


def normalize_mols(inputs: dict[str: Union[bool, list[Chem.Mol]]]) -> list[Chem.Mol]:
    """
    Take the input molecules and do the normalization steps as requested.
    """
    normalized_mols = []
    uncharger = rdMolStandardize.Uncharger()
    metal_disconnector = rdMolStandardize.MetalDisconnector()
    for mol in inputs['mols']:
        if mol is None:
            normalized_mols.append(None)
            continue

        norm_mol = Chem.Mol(mol)
        norm_mol = rdmolops.RemoveHs(norm_mol)
        if inputs['normalize']:
            norm_mol = rdMolStandardize.Normalize(norm_mol)
        if inputs['neutralize']:
            # we have been promised a helper function for this!
            norm_mol = uncharger.uncharge(norm_mol)
        if inputs['reionize']:
            norm_mol = rdMolStandardize.Reionize(norm_mol)
        if inputs['parent_fragment']:
            norm_mol = rdMolStandardize.FragmentParent(norm_mol)
        if inputs['metal_disconnect']:
            norm_mol = metal_disconnector.Disconnect(norm_mol)
        if inputs['canonical_tautomer']:
            norm_mol = rdMolStandardize.CanonicalTautomer(norm_mol)

        # Do a sanitize, as some of the steps don't do this, it seems.
        sanitize_mol(norm_mol)
        rdmolops.AssignStereochemistry(norm_mol)
        normalized_mols.append(norm_mol)

    return normalized_mols


def extract_input_values(request: DataFunctionRequest) -> tuple[dict[str: Union[bool, list[Chem.Mol]]],
                                                                str, DataType]:
    inputs = {}
    column_id = string_input_field(request, 'structureColumn')
    input_column = request.inputColumns[column_id]
    inputs['input_column_name'] = input_column.name
    mols = column_to_molecules(input_column)
    inputs['mols'] = mols

    for ids in [('normalizeID', 'normalize'),
                ('neutralizeID', 'neutralize'),
                ('reionizeID', 'reionize'),
                ('parentID', 'parent_fragment'),
                ('metalID', 'metal_disconnect'),
                ('tautomerID', 'canonical_tautomer')]:
        inputs[ids[1]] = boolean_input_field(request, ids[0])

    return inputs


def execute(request: DataFunctionRequest) -> DataFunctionResponse:
    inputs = extract_input_values(request)
    normalized_mols = normalize_mols(inputs)
    output_column = molecules_to_column(normalized_mols,
                                        f'Normalized {inputs["input_column_name"]}',
                                        DataType.BINARY)
    response = DataFunctionResponse(outputColumns=[output_column])
    return response
