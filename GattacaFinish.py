import sys
import os
from collections import OrderedDict
import argparse
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3
import JavaCodeGen
from gattaca_postplots import Plots
import glob
from scipy.stats import binom, gamma, poisson
import pickle
import time
from gattaca_classes import GeneLocs, FastaRecords
import multiprocessing as mp
import numpy as np
from numba import jit
from tqdm import tqdm
import psutil

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()

    requiredNamed = parser.add_argument_group("Required Arguments")
    inFile = parser.add_argument("-i", "--file", dest="inFile", default="./tests/stemCA.3.0/", help="Output file or directory from simulations.", required=True) #TODO change to true
    output = parser.add_argument("-o", "--output", dest="output", default="./tests/stemCA.3.0/gattaca/", help="Destination for directory to be created.", required=True)
    geneList = parser.add_argument("-f", "--geneList", dest="geneList", default="./tests/TestGenes.txt", help="List of genes, one per line to use for GattacaExample.Gattaca.", required=True)
    genome = parser.add_argument("-g", "--genome", dest="genome", default="GRCh37.75", help="Reference genome to use.")
    depth = parser.add_argument("-d", "--depth", dest="depth", default=500, help="Sequencing depth. Default 500.")
    depthDist = parser.add_argument("-k", "--depthFile", dest="depthFile", default=None, help="If a distribution for depth is desired use a file with one depth per variant per line.")
    cutoff = parser.add_argument("-c", "--cutoff", dest="cutoff", default=0.001, help="Minimum variant allele frequency (higher depth=lower cutoff recommended) This cutoff is for the Raw VAF, uncorrected. Default is a generous 0.001.")
    name = parser.add_argument("-n", "--name", dest="name", default="gattacaRun", help="Name of run. Default=gattacaRun.")
    tmpts = parser.add_argument('-t', '--timepoints', dest="tmpts", nargs='+', help='List of timepoints. Leave blank if you only want the last timepoint.', required=False)
    plots = parser.add_argument("-p", "--getPlots", dest="plots", default=True, action='store_false', help="Specifies whether to create plots.")
    vcfs = parser.add_argument("-vcf", "--getvcfs", dest="vcfs", default=False, action='store_true', help="Whether or not to construct VCF files from mutations.")
    verbose = parser.add_argument("-v", "--verbose", dest="verbose", default=False, action='store_true', help="Whether verbose output is turned on (development/debugging).")

    Options = parser.parse_args()  # main user args

    Options.output=os.path.expanduser(Options.output)
    Options.inFile=os.path.expanduser(Options.inFile)

    if Options.output[len(Options.output)-1]!="/":
        Options.output+="/"
    if Options.inFile[len(Options.inFile)-1]!="/" and Options.inFile[len(Options.inFile)-3:]!="csv":
        Options.inFile+="/"

    # try:
    #     Options.mutRate = float(Options.mutRate)
    # except:
    #     sys.exit("Improper mutation rate argument. Ensure this can be parsed as a decimal.")
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

class Data:

    def __init__(self, Options):
        # Raw Information
        self.data = Options.inFile
        print("Processing %s"%(self.data))
        self.tmpts = None
        self.cloneIDs = None
        self.parentIDs = None
        self.AllMuts = {}
        self.parsedData = OrderedDict() # Holds clone information with full mutational genome
        if Options.verbose:
            print("Data has been parsed.")
        self._readFile()
        if Options.verbose:
            print("File has been read. Constructing lineages.")
        self.lineages = {}
        start_time = time.process_time()
        self._smartLineages()
        if Options.verbose:
            print("Lineages have been constructed.")
            print("Finished building lineages: %.8f seconds" % (time.process_time() - start_time))
        self._getGenomes()
        if Options.verbose:
            print("Genomes have been built.")
            print("Finished building genomes: %.8f seconds" % (time.process_time() - start_time))

        self._getMutationsAtTimepoint()

        # Stats
        self.contexts = self._getOrder()

        # self.mutsAtTimepoint = []
        self.clonesWithMutsAtTimepoint = []
        self.popsWithMutsAtTimepoint = []
        self.popSize = [] # Built in self._getMutContexts(); list of values at each timepoint
        self.mutContexts = [] # list of context lists at each timepoint
        self.VAFs = []

        self._getMutContexts()
        self._getRawVAFs(Options)

    def _readFile(self):
        with open(self.data, "r") as theFile:
            lines = theFile.readlines()

        lines = [line.replace('\n','').split(',') for line in lines]
        idxTmptsStart = [i for i, k in enumerate(lines[0]) if k.startswith("ParentID")][0]+1

        self.tmpts = [int(float(i)) for i in lines[0][idxTmptsStart:len(lines[0])]]

        self.cloneIDs = []
        self.parentIDs = []
        for i, vals in enumerate(lines[1:]):
            self.cloneIDs.append(vals[4])
            self.parentIDs.append(vals[5])
            self.parsedData.update({ int(vals[4]) : {'CloneID': vals[4], 'hsv':(vals[1],vals[2],vals[3]), 'ParentID': vals[5], 'Pops': [i for i in vals[idxTmptsStart:len(vals)]], 'Genome': vals[0], 'NumMuts': 0 }  })


    def _getClonesPrivateGenome(self, cloneID):
        return(self.parsedData[int(cloneID)]['Genome'])

    def _smartLineages(self):
        '''
        Constructs lineages without having to iterate through the entire tree as a lot of clones will have info in childrens tree
        :return:
        '''
        clone_dict = dict(zip(self.cloneIDs, self.parentIDs))

        for i in tqdm(range(len(self.cloneIDs)), desc="Building lineages"):
            c = self.cloneIDs[i]
            p = clone_dict[c]
            lineage = [c, p]
            while p!='0':
                p = clone_dict[p] # get parents parent
                lineage.append(p)
            self.parsedData[int(c)].update({'Lineage':lineage})

        clone_dict = None

    def _getGenomes(self):
        '''
        Constructs the full genome of this clone based off the lineage information
        :return:
        '''
        clone_genome_dict = dict(zip(self.cloneIDs, [self.parsedData[int(c)]['Genome'] for c in self.cloneIDs]))

        processed = []
        for i in tqdm(range(len(self.cloneIDs)), desc="Building genomes; Mem=%s"%(psutil.virtual_memory().percent)):
            c = self.cloneIDs[i]
            genome = ''
            l = self.parsedData[int(c)]['Lineage']
            for ind_genome in l:
                genome += clone_genome_dict[ind_genome]
            self.parsedData[int(c)]['Genome']=genome

        for val in list(set(clone_genome_dict.values())):
            for v in val.split(';'):
                self.AllMuts.update({v:{'Pop':dict(zip(self.tmpts,[0 for i in self.tmpts]))}})

    def _getMutationsAtTimepoint(self):
        '''
        Gets the mutations that are at every timepoint and their pop at that timepoint
        :return:
        '''
        for i in tqdm(range(len(self.cloneIDs)), desc="Building time data; Mem=%s"%(psutil.virtual_memory().percent)):
            c = int(self.cloneIDs[i])
            info = self.parsedData[c]
            alive = [self.tmpts[int(i)] for i, val in enumerate(info['Pops']) if int(val)>0]
            pops = [int(val) for i, val in enumerate(info['Pops']) if int(val)>0]
            for mut in info['Genome'].split(';'):
                for i, val in enumerate(info['Pops']):
                    if int(val)>0:
                        self.AllMuts[mut]['Pop'][self.tmpts[int(i)]] += int(val)

    def _getOrder(self):
        bases = ['A','C','G','T']
        muts = ['[C>A]','[C>G]','[C>T]','[T>A]','[T>C]','[T>G]']
        vals = []
        for pos2 in range(0, 6):
            for pos1 in range(0,4):
                for pos3 in range(0,4):
                        vals.append(bases[pos1]+muts[pos2]+bases[pos3])
        return(vals)

    def _getMutContexts(self):
        for i, t in enumerate(self.tmpts):
            self._totalPop(i)

            out = OrderedDict.fromkeys(self.contexts)
            for k in out:
                out[k] = 0
            self.mutContexts.append(out)

        # Place counts based on all muts
        for mut in self.AllMuts:
            for i, tmpt in enumerate(self.tmpts):
                if self.AllMuts[mut]['Pop'][tmpt]!=0 and mut!='':
                    self.mutContexts[i][mut.split('.')[3]]+=1

    def _totalPop(self, timeIdx):
        val = 0
        for c in self.cloneIDs:
            val += int(self.parsedData[int(c)]['Pops'][timeIdx])
        self.popSize.append(val)

    def _getRawVAFs(self, Options):
        '''
        Goal, obtain the raw variant allele frequencies based on population numbers and individual numbers of cells with that mutation.
        :param Options:
        :return:
        '''
        muts = list(self.AllMuts.keys())
        for i in tqdm(range(len(muts)), desc="Calculating variant allele frequencies; Mem=%s" % (psutil.virtual_memory().percent)):
            thisMut = muts[i]
            self.AllMuts[thisMut].update({'Seq':{}}) # Places a timepoint key
            for idx, t in enumerate(self.AllMuts[thisMut]['Pop']):
                if self.AllMuts[thisMut]['Pop'][t]!=0:
                    n = self.popSize[idx]
                    vaf = self.AllMuts[thisMut]['Pop'][t]/(2.*float(n))
                    assert vaf<=1.0, "VAF Greater than 1.0 error: %s"%(vaf)
                    reads,depth,corrected_vaf = self._getCorrectedVAF(Options, vaf)
                    self.AllMuts[thisMut]['Seq'].update({t:{'rawVAF': vaf, 'correctedVAF': corrected_vaf, 'reads':reads,'depth':depth}})

    def _getCorrectedVAF(self, Options, vaf):
        '''
        Gets a corrected VAF based on the rawVAF using either a poisson distribution w/
        binomial or a gamma fit dist if options are selected
        :param Options: Command arguments
        :return: Corrected VAF, Reads, Depth
        '''
        if Options.depthFile != None:
            pass # TODO need to add GAMMA distribution fits from a file
        else:
            return(_getFastCorrectedVAF(Options.depth, vaf))
            # reads = binom.rvs(Options.depth, vaf)
            # depth = float(poisson.rvs(Options.depth))
            # cvaf = reads/depth
            # assert cvaf <= 1.0, "VAF Greater than 1.0 error: %s" % (cvaf)
            # return([reads,depth,cvaf])

@jit(nopython=True)
def _getFastCorrectedVAF(depth, vaf):
    d = float(np.random.poisson(depth))
    reads = np.random.binomial(depth, vaf)
    return(reads, d, reads/d)


def BuildOutputTable(Options, data, snpeff, refGenome):
    geneLocs = GeneLocs(Options, snpeff)
    fastaRec = FastaRecords(geneLocs, Options, snpeff, refGenome)
    revDict = {'t': 'a', 'a': 't', 'g': 'c', 'c': 'g'}

    outline = []
    for i, sim in enumerate(data):
        if Options.tmpts==None: # Means we are getting every tmpt
            for m in sim.AllMuts:
                mut = sim.AllMuts[m]
                for t in mut['Pop']:
                    if mut['Pop'][t]>0 and m!='':
                        vafData = mut['Seq'][t]
                        if mut['Seq'][t]['rawVAF']>Options.cutoff:
                            # Sample    Tmpt    Chrom   Pos Ref Alt Gene    Induction_t Trinuc TrinucTrue  RawVAF  CorrVAF Reads Depth
                            mutInfo = m.split('.')
                            ref = mutInfo[3].split('[')[1].split('>')[0]
                            alt = mutInfo[3].split('[')[1].split('>')[1].split(']')[0]
                            val = [os.path.basename(sim.data).replace('.csv',''), t, mutInfo[4], mutInfo[5], ref, alt, mutInfo[2], mutInfo[0], mutInfo[3][0] + alt + mutInfo[3][6], mutInfo[3][0] + alt + mutInfo[3][6], vafData['rawVAF'], vafData['correctedVAF'], vafData['reads'], vafData['depth']]
                            outline.append(val)
        else: # TODO build so that you get only specificed tmpts
            pass

    outFinal = []
    for mutOut in outline:
        for rec in fastaRec.geneSeq:
            if rec[3]==mutOut[6]:
                start = int(rec[1])
                end = int(rec[2])

                if mutOut[4].lower()==rec[0][int(mutOut[3])-start-1]:
                    outFinal.append(mutOut)
                else:
                    mutOut[4]=rec[0][int(mutOut[3])-start-1].upper()
                    mutOut[5]=revDict[mutOut[5].lower()].upper()
                    mutOut[9]=revDict[mutOut[9][0].lower()].upper()+revDict[mutOut[9][1].lower()].upper()+revDict[mutOut[9][2].lower()].upper()
                    outFinal.append(mutOut)

    WriteVCF(Options, outFinal)
    vcfData = AnnotateVCF(Options, snpeff, Options.output + "tmp.vcf")

    anotherList = []
    for mutOut in outFinal:
        key = mutOut[2] + "\t" + mutOut[3] + "\t" + mutOut[4] + "\t" + mutOut[5]
        effectField = vcfData[key].split('\t')[7].split(";")[2].split("=")[1].split("|")
        effect = effectField[1]
        mutOut.append(effect)
        anotherList.append(mutOut)

    with open(Options.output + "mutations.txt", "w") as outfile:
        outfile.write("Sample\tTmpt\tChrom\tPos\tRef\tAlt\tGene\tInduction_t\tTrinuc\tTrinucTrue\tRawVAF\tCorrVAF\tReads\tDepth\tEff\n")
        for item in anotherList:
            item = [str(i) for i in item]
            outfile.write("\t".join(item) + "\n")

    print("Mutation table output to: %s"%(Options.output + "mutations.txt"))

def BuildOutputCloneData(Options, data):
    outline = []
    for i, sim in enumerate(data):
        for clone in sim.parsedData:
            line = [os.path.basename(sim.data).replace('.csv',''), sim.parsedData[clone]['CloneID'], sim.parsedData[clone]['ParentID'], len(sim.parsedData[clone]['Genome'].split(';')), '\t'.join(sim.parsedData[clone]['Pops'])]
            outline.append('\t'.join([str(item) for item in line]))

    with open(Options.output + "mutations.clonal.txt", "w") as outCounts:
        outCounts.write('\t'.join(["Sample","CloneID", "ParentID", "NumMuts",'\t'.join([str(i) for i, v in enumerate('\t'.join(sim.parsedData[0]['Pops']))])]) + "\n")
        outCounts.write("\n".join(outline))
    print("Clone data table output to %s"%(Options.output + "mutations.clonal.txt"))

def AnnotateVCF(Options, snpeff, inFile):
    print("Annotating mutations...")
    cmd = ["java", "-jar", snpeff['snpeff'], "ann", "-canon", Options.genome, inFile]
    locInfo = os.popen(' '.join(cmd)).read()
    with open(inFile.replace('.vcf','.ann.vcf'), 'w') as outfile:
        outfile.write(locInfo)
    print("Annotation Complete...file written to %s"%(inFile.replace('.vcf','.ann.vcf')))
    os.remove(inFile)

    out = {}
    for line in locInfo.split('\n'):
        if line.startswith("#")==False and line!='':
            k = line.split('\t')
            key = k[0] + "\t" + k[1] + "\t" + k[3] + "\t" + k[4]
            out.update({ key : line})
    return(out)

def WriteVCF(Options, mutOut):
    header = "##fileformat=VCFv4.2\n##fileDate=%s\n##source=GattacaExample.Gattaca\n##reference=%s\n##contig=None\n##phasing=partial"%(time.asctime( time.localtime(time.time()) ),Options.genome)
    header2 = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
    header3 = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Variant Read Depth\">"
    header4 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tALLSAMPLES"
    seen = []
    with open(Options.output + "tmp.vcf", 'w') as outFile:
        outFile.write('\n'.join([header,header2,header3,header4]) + "\n")
        for m in mutOut:
            out = [m[2],m[3],".",m[4],m[5],".","PASS"]
            if out not in seen:
                infoField = "DP=%s;AF=%s\tGT:DP:VR\t0/1:%s:%s"%(int(m[13]),round(m[11],5),int(m[13]),m[12])
                outFile.write('\t'.join(out) + "\t" + infoField + "\n")
                seen.append(out)

if __name__=="__main__":
    # Pull information about usr system files
    localpath = os.path.abspath(__file__).replace('GattacaFinish.py', '')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Options = Parser()
    Config.read(localpath + "usr_paths.ini")
    snpeff = ConfigSectionMap(Config.sections()[0], Config)  # get annovar script paths
    refGenome = ConfigSectionMap(Config.sections()[1], Config)  # get reference genome

    if os.path.isdir(os.path.expanduser(Options.output))==False:
        os.mkdir(os.path.expanduser(Options.output))
        if Options.verbose:
            print("Output directory created.")
    else:
        if Options.verbose:
            print("Output drictory found.")

    data = []
    if os.path.isfile(os.path.expanduser(Options.inFile)): # Individual files
        if Options.verbose:
            print("Creating data class.")
        data.append(Data(Options))
    elif os.path.isdir(os.path.expanduser(Options.inFile)): # Multiple files
        allFiles = glob.glob(Options.inFile + "*.csv")
        Options2 = Options
        if Options.verbose:
            print("Directory found, creating data classes.")
        for i in allFiles:
            Options2.inFile = i
            data.append(Data(Options2))
    else:
        sys.exit("No input file(s) found.")
    #
    # # # TODO build endpoint table or table with all mutations.
    # pickle.dump(data, open('postProcess.p', 'wb')) # TODO Delete when done with dev
    # sys.exit()
    # data = pickle.load(open('postProcess.p', 'rb')) # TODO Delete when done with dev
    BuildOutputCloneData(Options, data)
    BuildOutputTable(Options, data, snpeff, refGenome)
    sys.exit("Here")
    if Options.plots:
        build = Plots(Options,data)
        # build.mutContextMovie(Options, data)
        # build.mutVAF(Options, data)
        build.oneOverF(Options, data)



