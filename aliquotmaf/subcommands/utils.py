"""Ulitity functions"""
import gzip
import re
import json

def get_open_function(fil):
    """
    Returns the appropriate open function based on gz ending.
    """
    if fil.endswith('.gz'):
        return gzip.open
    return open

def load_json(fil):
    """
    Loads a json file.
    """
    dat = None
    with open(fil, 'rt') as fh:
        dat = json.load(fh)
    return dat

def assert_sample_in_header(vcf_object, sample, can_fail=False):
    """
    Asserts that a given sample is in the VCF header and returns the index.

    :param vcf_object: ``~pysam.VariantFile`` instance
    :param sample: the sample's name
    :return: sample column index
    """
    idx = None
    if not can_fail:
        assert sample in vcf_object.header.samples, \
            'Unable to find sample {0} in VCF header!'.format(sample)

    # Get the index of the samples
    try:
        slist = list(vcf_object.header.samples)
        idx = slist.index(sample) + 9
    except ValueError:
        pass
    return idx

def extract_annotation_from_header(vcf_object, vep_key='CSQ'):
    """
    Extract the VEP annotation columns from the VCF header.

    :param vcf_object: ``~pysam.VariantFile`` instance
    :return: annotation column list
    :return: the key in the INFO data from the VCF
    """
    # annotation cols
    ann_cols_format = []

    # Loop over header records
    for record in vcf_object.header.records:
        if record.type == 'INFO':
            iname = record.get('ID')
            if iname and str(iname) == vep_key:
                vep_key = iname
                anno_line = re.search('Format: (\S+)"$', record['Description'])
                raw_ann_cols_format = anno_line.group(1).split('|')
                for ann in raw_ann_cols_format:
                    if ann.startswith('ExAC'): ann_cols_format.append(fix_exac(ann))
                    else: ann_cols_format.append(ann)

    assert vep_key is not None, 'Malformed header, unable to find VEP annotation key'
    assert ann_cols_format, 'Malformed header, unable to find VEP annotation column info'
    return ann_cols_format, vep_key

def fix_exac(ann):
    """
    Convert the to the "OLD" version of ExAC headers from the old plugin.
    before VEP included it built in.

    :param ann: the annotation key
    :return: the converted annotation key
    """
    lookup = {"ExAC_MAF" : "ExAC_AF",
        "ExAC_AFR_MAF" : "ExAC_AF_AFR",
        "ExAC_AMR_MAF" : "ExAC_AF_AMR",
        "ExAC_EAS_MAF" : "ExAC_AF_EAS",
        "ExAC_FIN_MAF" : "ExAC_AF_FIN",
        "ExAC_NFE_MAF" : "ExAC_AF_NFE",
        "ExAC_OTH_MAF" : "ExAC_AF_OTH",
        "ExAC_SAS_MAF" : "ExAC_AF_SAS",
        "ExAC_Adj_MAF" : "ExAC_AF_Adj"}
    return lookup[ann]

def load_enst(fpath):
    """
    Loads the custom transcript overrides file if the user provided it

    :param fpath: the custom override file

    :returns set: a set of ENST identifiers
    """
    lst = []
    with open(fpath, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            cols = line.rstrip('\r\n').split('\t')
            lst.append(cols[0])
    return set(lst)
