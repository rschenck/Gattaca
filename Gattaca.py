'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''

import sys
import os
import argparse
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3
import JavaCodeGen
import pickle

from gattaca_classes import GeneLocs, FastaRecords, Contexts
from gattaca_muts import Mutations
from gattaca_plots import MutPropPlots
from JavaCodeGen import JavaCode


def Parser():
    # get user variables
    parser = argparse.ArgumentParser()

    requiredNamed = parser.add_argument_group("Required Arguments")
    geneList = parser.add_argument("-f", "--geneList", dest="geneList", default="./tests/TestGenes.txt", help="List of genes, one per line to use for GattacaExample.Gattaca.", required=True)
    output = parser.add_argument("-o", "--output", dest="outputLoc", default="./tests/GattacaEx/src/", help="Output directory for classes to be used with HAL.", required=True)
    genome = parser.add_argument("-g", "--genome", dest="genome", default="GRCh37.75", help="Reference genome to use.")

    outputNamed = parser.add_argument_group("Output options")
    outputNamed.add_argument("-b", "--getBams", dest="getBams", default=False, action='store_true', help="Specifies whether to create a BAM file for use post simulation.")

    intputNamed = parser.add_argument_group("Input")
    intputNamed.add_argument("-c", "--contextFile", dest="contextFile", default="",
                        help="Trinucleotide context fractions, properly formatted. See github repo.")

    mutNamed = parser.add_argument_group('Mutation Information')
    mutNamed.add_argument("-u", "--mutRate", dest="mutRate", default="3.2E-9", help="Mutation rate. Default 3.2E-9.")

    Options = parser.parse_args()  # main user args

    if Options.output[len(Options.output)-1]!="/":
        Options.output+="/"

    try:
        Options.mutRate = float(Options.mutRate)
    except:
        sys.exit("Improper mutation rate argument. Ensure this can be parsed as a decimal.")
    return(Options)

def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

if __name__=="__main__":
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('Gattaca.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Options = Parser()
    Config.read(localpath + "usr_paths.ini")
    snpeff = ConfigSectionMap(Config.sections()[0], Config)  # get annovar script paths
    ref = ConfigSectionMap(Config.sections()[1], Config)  # get reference genome

    '''
    Step 1: Read in the target names. Get the gene information (length) location. etc.
    '''
    geneLocs = GeneLocs(Options, snpeff)
    pickle.dump(geneLocs, open('file.p', 'wb')) # TODO Delete when done with dev
    geneLocs = pickle.load(open('file.p', 'rb')) # TODO Delete when done with dev

    '''
    Step 2: Get the sequences into a FASTA file
    '''
    fastaRec = FastaRecords(geneLocs, Options, snpeff, ref)

    for g in fastaRec.geneSeq:
        if len(g[0])==0:
            print("No gene found: " + g[3])

    '''
    Step 3: Get the contexts and probabilities
    '''
    contexts = Contexts(Options, geneLocs, fastaRec)

    '''
    Step 4: Get mutation information for genes
    '''
    mutRates = Mutations(geneLocs.data, Options.mutRate)

    '''
    Step Optional: Get some plots
    '''
    # plots = MutPropPlots(contexts.context, contexts.triNucsProb, contexts.triNucsGenePos, contexts.context32Order, mutRates, geneLocs)# TODO uncomment when done with dev

    '''
    Last Step:
    '''
    JavaCode([gene for gene in geneLocs.data], contexts.context, mutRates.geneMuts, contexts.triNucsProb, contexts.triNucsGenePos, Options.outputLoc)