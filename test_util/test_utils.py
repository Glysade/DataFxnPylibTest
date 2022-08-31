from typing import Dict

def service_class_to_name_dict() -> Dict[str, str]:
    return {
        'BlastWebService': 'Blast (web)',
        'BlastTableService': 'Blast (table)',
        'BlastLocalSearchService': 'Blast (local) column query',
        'BlastLocalSequenceSearchService': 'Blast (local) text query',
        'SequenceAlignService': 'Multiple Sequence Alignment',
        'Omega2ResourceService': 'Omega2 (resource)',
        'Omega2Service': 'Omega2',
        'RocsResourceSearchService': 'ROCS Database Search',
        'RocsTableSearchService': 'ROCS Table Search',
        'ExtractFeatureService': 'Extract sequence feature',
        'SubsequenceService': 'Extract sequence range',
        'RgroupDecompositionService': 'R Group Decomposition',
        'MmpdbService': 'MMPDB Database Search',
        'PairwiseAlignmentService': 'Pairwise Sequence Alignment',
        'SequenceRetrievalService': 'Sequence Retriever',
        'AntibodyNumberingService': 'Antibody Numbering',
        'ConsensusSequenceService': 'Consensus sequence'
    }


def ruse_service_name(service_class: str) -> str:
    """Maps service class to current ruse service name"""

    service_class_to_name = service_class_to_name_dict()

    if service_class not in service_class_to_name:
        raise ValueError("no ruse task defined for class {}".format(service_class))
    return service_class_to_name[service_class]


