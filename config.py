# place to download and get, UCSC metadata.
KNOWN_GENES = 'meta_data/clustersWithNames.tsv'


# the data itself - directory of mat files (matlab) with convention REPLICATE_CONDITION_TIME.mat
# where:
#   REPLICATE is some unique value (such as 384)
#   CONDITION is some short description of condition (such as control)
#   TIME is duration we collect the cells (relative to time we washed out the DRB)
TRANSCRIPTION_DATA_DIR = 'mat/'

# where to write the output to
HMM_RESULT_DIR = 'hmm_raw_res/'

ALL_REPLICATES = []  # replicates (OPTIONAL - if you sue the script without specifying replicate)

# ======= OTHER CONFIGURATIONS
MIN_GENE_LENGTH = 35000  # minimum length of genes to take into account. (short genes will be ignored)
jump = 50  # jumps of bins in the matlab array/resolution
min_sequence_length = 500 / jump  # minimum sequence length after removing exons. shorter transcripts will be skipped
