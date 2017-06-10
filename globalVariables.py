# globalVariables.py

MIN_OVERLAP = ''
TIP_SIZE = ''
CONFIDENCE_BASE = ''
ERROR_CORRECTION_THRESHOLD = ''
READ_LEN = ''
FULL_LENGTH = ''
PATH_FINDING_PARAMETER = ''
NEED_DEFLANK = ''

# Output variables
graph_path = ''
plot_dir = ''
rj_dir = ''
full_genes_path = ''
fragments_path = ''
filename = ''

# Global dictionaries
read_db = {}
r_pos = {}
cm_pos = {}
read_db_original = {}
read_position_db = {}

def init():
# Input variables
    global MIN_OVERLAP
    global TIP_SIZE
    global CONFIDENCE_BASE
    global ERROR_CORRECTION_THRESHOLD
    global READ_LEN
    global FULL_LENGTH
    global PATH_FINDING_PARAMETER
    global NEED_DEFLANK

    # Output variables
    global graph_path
    global plot_dir
    global rj_dir
    global full_genes_path
    global fragments_path
    global filename
    global output_dir

    # Global dictionaries
    global read_db
    global r_pos
    global cm_pos
    global read_db_original
    global read_position_db
