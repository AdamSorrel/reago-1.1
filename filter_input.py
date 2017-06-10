import sys
import os
import argparse

###################################################################################
# User command line input
###################################################################################
# This help message can be called with the -h or --help option (e.g. filter_input.py -h)

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='# This program filters bacterial sequences from an input file using the cmsearch algorithm.')

parser.add_argument('-r1', '--read1', dest='R1', help='Input forward (1) read fasta file.')
parser.add_argument('-r2', '--read2', dest='R2', help='Input reverse (2) read fasta file.')
parser.add_argument('-o', '--output_dir', dest='OUT', default = 'Output', help='Output directory. [Default: Output]')
#parser.add_argument('-cmd', '--cm_directory', dest='CMDIR', default='cm' ,help='Directory with the cm files. [Default: cm]')
parser.add_argument('-cm', '--cm', dest='CM', default='cm/bacteria.cm' ,help='Relative path of the cm database to be used. [Default: cm/bacteria.cm]')
parser.add_argument('-cpu', '--cpu_number', dest='CPU', default='1' ,help='Number of cpu for parallelisation. [Default: 1]')

args = parser.parse_args()

###################################################################################
# Function definitions
###################################################################################

def get_rc(read):   # This function retrieves a reverse complement of the sequence
    read = read.upper()
    read = list(read)

    for i in range(len(read)):
        if read[i] not in "ACTGU":
            read[i] = "N"

    read = "".join(read)

    # DNA
    if read.find("T") >= 0:
        alphabet = {'A': 'T',
                    'G': 'C',
                    'C': 'G',
                    'T': 'A',
                    'N': 'N'}
    # RNA
    else:
        alphabet = {'A': 'U',
                    'G': 'C',
                    'C': 'G',
                    'U': 'A',
                    'N': 'N'}
    rev_com = ""
    for char in read:
        rev_com += alphabet[char]

    rev_com = rev_com[::-1]

    return rev_com



def get_fa(fn):
    seq_d = {}  # Sequence dictionary
    first = True    # Switch to skip the first sequence
    f = open(fn)
    for line in f:
        if line[0] == ">":
            if first != True:
                seq_d[seq_id] = seq.upper() # Saving string as all upper-case
            else:   # This is to skip the first sequence
                first = False

            seq_id = line.split()[0][1:]    # Splits the sequece ID and retreaves the header
            seq = ""
        else:
            seq += line[:-1]

    seq_d[seq_id] = seq.upper()
    f.close()

    return seq_d, seq_id[-1]

# Constructs a database from the cmsearch output
def process_single_file(out_fn):
    d = {}  # Database
    start = False   # This is to skip the first line
    #f = open(out_fn)
    with open(out_fn) as fIn: # Opening file in a proper python way
        for line in fIn:
            # This is in fact obsolete, since the lines are selected based on the '>>', which only appears in the verbose output.
            if line[:15] == "Hit alignments:":  # Scrolling through the entire beginning of the file all the way to the verbose output.
                start = True
            if start == False:
                continue

            if line[:2] == ">>":    # Saving sequence code
                read_id = line.split()[1]
                base = read_id[:-1]
            elif "!" in line:   # Retrieving parameters
                data = line.split()
                # m_st - template start position
                # m_ed - template end position
                # s_st - query start position
                # s_ed - query end position
                # strand - strand (+ or -)
                m_st, m_ed, s_st, s_ed, strand = data[6], data[7], data[9], data[10], data[11]
                if int(s_st) > int(s_ed):   # Making sure the query is always from smaller to larger number
                    s_st, s_ed = data[10], data[9]
                d[base] = [m_st, m_ed, s_st, s_ed, strand]  # Saving everything in a database
    #f.close()  # Unnecessary with files opened by the method "with"
    return d



# This function only performs fun process_single_file twice.
def process(out_fn_1, out_fn_2):
    d_1 = process_single_file(out_fn_1)
    d_2 = process_single_file(out_fn_2)
    return d_1, d_2

###################################################################################
# Retreaving output
###################################################################################
# All arguments are stored in a format 'args.dest', where dest is a destination defined in the argparse command
# e.g. args.OUT is a user defined output

fn_1 = args.R1  # User input forward read 1
fn_2 = args.R2  # User input reverse read 2
out_dir = args.OUT  # User defined output directory
cm = args.CM # User defined cm database
cpu = args.CPU # User defined number of CPU.

#args = sys.argv # Getting user input files
#try:
#    fn_1, fn_2, out_dir, cm_dir, cm, cpu = args[1:]

# Adding slash to the output directory
if out_dir[:-1] != "/":
    out_dir += "/"

# Adding slash to the output directory
#if cm_dir[:-1] != "/":
#    cm_dir += "/"

# Testing if the forward read file exists
if not os.path.exists(fn_1):
    print("Error:", fn_1, "doesn't exist.")
    sys.exit(1)

# Testing if the reverse read file exists
if not os.path.exists(fn_2):
    print("Error:", fn_2, "doesn't exist.")
    sys.exit(1)

# If an output folder doesn't exist, create it
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Generating names for the output from the cm database intpu
cmName = cm.split("/")[-1]  # Splitting string by / to bet rid of subfolders and getting the last element of the
cmName = cmName.split(".")[0]  # Splitting the file name to get rid of file type suffix (.cm in this case) and getting the 1st element (name of the database)


###################################################################################
# cmsearch run
###################################################################################

# Message to keep user entertained
print("Indentifying 16S reads")

# Run of the cmsearch for forward read (1)
print("Forward read.")
os.system("cmsearch --cpu " + cpu + " " + cm + " " + fn_1 + " > " + out_dir + cmName + "_1.out")
# Run of the cmsearch for reverse read (2)
print("Reverse read.")
os.system("cmsearch --cpu " + cpu + " " + cm + " " + fn_2 + " > " + out_dir + cmName + "_2.out")

###################################################################################
#
###################################################################################

# Reverse complements and constructs a database of all the forward reads
db_1, end_symbol_1 = get_fa(fn_1)
# Reverse complements and constructs a database of all the reverse reads
db_2, end_symbol_2 = get_fa(fn_2)


#if "b" in cm:
#    b_1, b_2 = process(out_dir + "b_1.out", out_dir + "b_2.out")
#else:
#    b_1, b_2 = {}, {}

#if "a" in cm:
#    a_1, a_2 = process(out_dir + "a_1.out", out_dir + "a_2.out")
#else:
#    a_1, a_2 = {}, {}

# Creating two databases of from the input files
db1, db2 = process(out_dir + cmName + "_1.out", out_dir + cmName + "_2.out")

# Intersection of both fwd and rws read databases
dbIntersection = set(db1.keys()).intersection(set(db2.keys()))

#b_base = set(b_1.keys()).intersection(set(b_2.keys()))

#a_base = set(a_1.keys()).intersection(set(a_2.keys()))


# Read count
read_cnt = 1
# Opening output file
fo = open(out_dir + "filtered.fasta", "w")
# ???
included_read = set([])

for entry in dbIntersection:
    # Retrieving header values from each of the databases
    # 1. template start, 2. template end, 3. query start, 4. query end, strand (+/-)
    m_st_1, m_ed_1, s_st_1, s_ed_1, strand_1 = db1[entry]
    m_st_2, m_ed_2, s_st_2, s_ed_2, strand_2 = db2[entry]

    # Joining header values in a string
    pos_str_1 = " ".join([m_st_1, m_ed_1, s_st_1, s_ed_1])
    pos_str_2 = " ".join([m_st_2, m_ed_2, s_st_2, s_ed_2])

    # Should both reads be of the same strand orientation, nothing else is done (why?)
    if strand_1 == strand_2:
        print(pos_str_1)
        print(pos_str_2)
        continue

    # Entry is a sequence code (e.g. "50434992.ACA.")
    # end_symbol_1 or 2 are just the last symbols from the header of input file
    # signifying whether it's a forward (1) or reverse (2) read.
    read_id_1 = entry + end_symbol_1
    read_id_2 = entry + end_symbol_2

    # Retrieving sequences from appropriate database
    seq_1 = db_1[read_id_1]
    seq_2 = db_2[read_id_2]

    # If either strand is negative, sequence is reverse transcribed
    if strand_1 == "-":
        seq_1 = get_rc(seq_1)
    if strand_2 == "-":
        seq_2 = get_rc(seq_2)

    # This theoretically makes sure that reads that has been included in Bacteria will not be included in Archaea.
    # This should not happen!
    included_read.add(read_id_1)
    included_read.add(read_id_2)

    # Writing the output to a file.
    fo.write(">" + str(read_cnt) + ".1 " + pos_str_1  + "\n")
    fo.write(seq_1 + "\n")
    fo.write(">" + str(read_cnt) + ".2 "  + pos_str_2  + "\n")
    fo.write(seq_2 + "\n")
    # Reads are labeled by their respective numbers in the output.
    read_cnt += 1

fo.close()
