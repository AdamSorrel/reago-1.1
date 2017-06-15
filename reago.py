import time
import networkx as nx
#import pygraphviz as pgv
import argparse
from argparse import RawTextHelpFormatter

from reagoFunctions import *
from parseInput import parseInput

# TODO : phase out global variables
# Imports global settins file (defining global variables)
import globalVariables as g

import multiprocessing as mp


#==============================================================================
# # Speed profiling %%
#==============================================================================


start = time.time()

#==============================================================================
# # This is for the floats to work properly
#==============================================================================


class fRange(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return ('{0}-{1}'.format(self.start, self.end))

#==============================================================================
# # User command line input
#==============================================================================

# This help message can be called with the -h or --help option (e.g. filter_input.py -h)

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='# Reago is an assembly program using the filter_input.py generated file and building a read assembly.\n \
Dependencies: \n \
\t 1. Readjoiner 1.2 \n \
\t 2. Infernal 1.1.1 \n',\
formatter_class=RawTextHelpFormatter)

# Optional arguments
parser.add_argument('-ol', '--overlap', dest='OVERLAP', type = float,  default=None, help='Read overlap. [Default: 0.7]', choices=[fRange(0.5, 1)])
parser.add_argument('-e', '--error', dest='ERROR_CORRECTION_THRESHOLD', type = float, default=None, help='Error correction threshold. [Default: 0.05]', choices=[fRange(0, float(1))])
parser.add_argument('-t', '--tip_size', dest='TIP_SIZE',  type = int, default=None, help='Tip size. [Default: 30]', choices=[fRange(1, float("inf"))])
parser.add_argument('-p', '--path_finding', dest='PATH_FINDING_PARAMETER',  type = int, default=None, help='Path finding parameter. [Default: 10]', choices=[fRange(2, 11)])
parser.add_argument('-f', '--full_length', dest='FULL_LENGTH',  type = int, default=None, help='Full length of the gene. [Default: 1350]', choices=[fRange(1, float("inf"))])
parser.add_argument('-c', '--confidence_base', dest='CONFIDENCE_BASE',  type = int, default=None, help='confidence_base. [Default: 10]', choices=[fRange(2, 11)])


# Required arguments
requiredNamed = parser.add_argument_group('Required arguments')
requiredNamed.add_argument('-i', '--input', dest='IN', default = None, help='Input fasta file generated by filter_input.py.', required=False)
requiredNamed.add_argument('-o', '--output', dest='OUT', default=None,help='Output folder.[Default: sample_out/]', required=False)
requiredNamed.add_argument('-l', '--read_length', dest='READ_LENGTH', type = int, help='Read length.', required=False, choices=[fRange(1, float("inf"))])

args = parser.parse_args()

#==============================================================================
# # Input parse
#==============================================================================

# Defining global variables
# TODO : phase global variables out and replace by JSON settings object
g.init()

# Function from reagoFunctions.py which assigns global variables based on user input
variables = parseInput(args)

print('It took', time.time()-start, 'seconds to parse the input.')

#==============================================================================
# # Main procedure starts here
#==============================================================================

start = time.time()


print (timestamp(), "REAGO (v1.10) started...")
print ("Input file:", variables.filename)

print (timestamp(), "Reading input file...")
# Retreaving sequence database (read_db), query position database (r_pos) and template position database (cm_pos)

# Creating a redis database object

rserv = redis.Redis('localhost')

get_fa(variables, rserv)

read_db_original = dict(variables.read_db)   # Saving database for future use in fun get_assemblie
initialize_read_pos(variables, d) # Sets read position database to 0. Saving database for future use in fun get_assemblie
combine_duplicated_reads(d) # Dereplicating a database and saving the derep. headers separated by a '|' character.
#Saving database for future use in fun get_assemblie

print('It took', time.time()-start, 'seconds to prepare databases.')
start = time.time()


print (timestamp(), "Initializing overlap graph...")
# Generating a DiGraph with a readjoiner using a file 'graph'.
G = create_graph_using_rj("graph", variables, dat)
# Passing the DiGraph to to networkx
# more info in https://networkx.github.io/documentation/development/reference/generated/networkx.algorithms.components.weakly_connected.weakly_connected_component_subgraphs.html#networkx.algorithms.components.weakly_connected.weakly_connected_component_subgraphs)
subgraphs = nx.weakly_connected_component_subgraphs(G)

# Retrieving scaffold candidates and full genes from the subgraphs
print (timestamp(), "Recovering 16S rRNAs...")
# Starting lists of results
full_genes = []
scaffold_candidates = []

print('It took', time.time()-start, 'second to prepare main graph.')
start = time.time()

quit()

#num = 0
# Looping through the output of networkx

# From python 3.4 onwards
#nCPU = os.cpu_count()

#nCPU = 4
#manager = mp.Manager()

#pool = mp.Pool(nCPU)
#pool.map(subgraph_treatment, subgraphs)
#jobs.append(p)
#p.start(p)

for subgraph in subgraphs:
    full_genes,scaffold_candidates = subgraph_treatment(subgraph)
    #p = multiprocessing.Process(target=worker)
    #jobs.append(p)
    #p.start()


print('It took', time.time()-start, 'second to run through subgraphs.')

print (timestamp(), "Scaffolding on short 16S rRNA segments...")
scaf, remaining = scaffold(scaffold_candidates)
full_genes += scaf

print (timestamp(), "Write to Files...")
write_gene(full_genes)
write_frag(remaining)


print (timestamp(), "Done.")
print ("- Number of 16S rRNAs:", len(full_genes))
print ("- Full genes:", g.full_genes_path)
print ("- Gene fragments:", g.fragments_path)
