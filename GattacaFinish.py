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
    cutoff = parser.add_argument("-c", "--cutoff", dest="cutoff", default=0.005, help="Minimum variant allele frequency (higher depth=lower cutoff recommended). Default is a generous 0.005.")
    name = parser.add_argument("-n", "--name", dest="name", default="gattacaRun", help="Name of run. Default=gattacaRun.")
    tmpts = parser.add_argument('-t', '--timepoints', dest="tmpts", nargs='+', help='List of timepoints. Leave blank if you only want the last timepoint.', required=False)
    plots = parser.add_argument("-p", "--getPlots", dest="plots", default=True, action='store_false', help="Specifies whether to create a BAM file for use post simulation.")
    vcfs = parser.add_argument("-vcf", "--getvcfs", dest="vcfs", default=False, action='store_true', help="Whether or not to construct VCF files from mutations.")
    verbose = parser.add_argument("-v", "--verbose", dest="verbose", default=True, action='store_false', help="Whether verbose output is turned on (development/debugging).")

    Options = parser.parse_args()  # main user args

    Options.output=os.path.expanduser(Options.output)
    Options.inFile=os.path.expanduser(Options.inFile)

    if Options.output[len(Options.output)-1]!="/":
        Options.output+="/"
    if Options.inFile[len(Options.inFile)-1]!="/":
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
        self.parsedData = OrderedDict() # Holds clone information with full mutational genome
        self._readFile()
        self.lineages = {}
        start_time = time.clock()
        self._buildLineages()
        self._buildFullTimepointGenome()
        self.contexts = self._getOrder()
        print("Finished building lineages: %.2f seconds"%(time.clock()-start_time))

        # Stats
        '''Provides information on a per timepoint snapshot {time 1: {clone info}, time 2: {clone info}, etc.}'''
        self.mutsAtTimepoint = []
        self.clonesWithMutsAtTimepoint = []
        self.popsWithMutsAtTimepoint = []
        self.popSize = []
        self.mutContexts = []
        self.VAFs = []
        start_time = time.clock()
        self._getMutContexts()
        print("Finished extracting mutation contexts: %.2f seconds"%(time.clock()-start_time))
        start_time = time.clock()
        self._getRawVAFs(Options)
        print("Finished getting VAFs: %.2f seconds"%(time.clock()-start_time))
        # Get Muts per clone per timepoint

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


    def _buildFullTimepointGenome(self):
        vals = [0 for i in range(0,99)]
        for i, c in enumerate(self.cloneIDs):
            # print("Building full genome for clone: " + c)
            genome = self.parsedData[int(c)]['Genome']
            assert int(self.lineages[c][0])==int(c), "Clones may not match Clone: %s, Ancestor 1st Clone: %s (these two should be equal)."%(c, ancestor[0])

            for j, ancestor in enumerate([int(item) for item in self.lineages[c]]):
                if ancestor!=int(c) and ancestor != 0:
                    genome+=self._getClonesPrivateGenome(ancestor)

            self.parsedData[int(c)]['Genome']=genome
            self.parsedData[int(c)]['NumMuts']=len(genome.split(';'))

    def _getClonesPrivateGenome(self, cloneID):
        return(self.parsedData[int(cloneID)]['Genome'])

    def _buildLineages(self):
        for i, clone in enumerate(self.cloneIDs):
            lineage = [clone]
            p = self.parentIDs[i] # starting parent
            lineage.append(p)
            while p!='0': # continue until parent is 0 and clone is 0
                # parent of every parent.
                p = self._getCloneIDofParent(p)
                lineage.append(p)
            self.lineages.update({clone:lineage})

    def _getCloneIDofParent(self, parent):
        for i, clone in enumerate(self.cloneIDs):
            if clone==parent:
                return(self.parentIDs[i])

    def _IDoNotKnowWhatWasHere(self, clone):
        for k, parent in enumerate(self.cloneIDs):
            if parent==clone:
                return(parent)

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
            self._getMutationsAtTmpt(i)

            out = OrderedDict.fromkeys(self.contexts)
            for k in out:
                out[k] = 0

            for mut in self.mutsAtTimepoint[i]:
                out[mut.split('.')[3]]+=1

            self.mutContexts.append(out)
        # print(self.mutContexts)

    def _getMutationsAtTmpt(self, timeIdx):
        muts = []
        for c in self.cloneIDs:
            if int(self.parsedData[int(c)]['Pops'][timeIdx]) != 0:
                m = self.parsedData[int(c)]['Genome'].split(';')
                for singleMut in m:
                    if singleMut != '':
                        muts.append(singleMut)
        self.mutsAtTimepoint.append(list(set(muts)))

    def _totalPop(self, timeIdx):
        val = 0
        for c in self.cloneIDs:
            val += int(self.parsedData[int(c)]['Pops'][timeIdx])
        self.popSize.append(val)

    def _getRawVAFs(self, Options):
        start_time = time.clock()
        cellsWithMut = OrderedDict()
        for i, muts in enumerate(self.mutsAtTimepoint):
            theMuts={}
            for mut in muts:
                theMuts.update({mut: {'cells': 0.0, 'rawVAF': 0.0, 'correctedVAF': 0.0, 'reads': 0, 'depth': 0}})
            cellsWithMut.update({i:theMuts})
        if(Options.verbose):
            print("Setup muts: %.2f seconds"%(time.clock()-start_time))

        start_time = time.clock()
        for clone in self.parsedData: # Clone
            g = self.parsedData[int(clone)]['Genome'].split(';')
            for mut in g:
                if mut != '':
                    for i, val in enumerate(self.parsedData[int(clone)]['Pops']):
                        if int(val)>0:
                            cellsWithMut[i][mut]['cells']=int(val)

        if (Options.verbose):
            print("Cells extracted: %.2f seconds" % (time.clock() - start_time))


        # for i, muts in enumerate(self.mutsAtTimepoint):
        #     theMuts={}
        #     for mut in muts:
        #         theMuts.update({mut: {'cells': 0.0, 'rawVAF': 0.0, 'correctedVAF': 0.0, 'reads': 0, 'depth': 0}})
        #         thePop = 0
        #         for c in self.clonesWithMutsAtTimepoint[i]:
        #             if mut in self.parsedData[int(c)]['Genome']:
        #                 thePop+=int(self.parsedData[int(c)]['Pops'][i])
        #         theMuts[mut]['cells']=thePop
        #
        #     cellsWithMut.update({i:theMuts})

        start_time = time.clock()
        for t, tmpt in enumerate(self.tmpts):
            for mut in cellsWithMut[t]:
                try:
                    cells = cellsWithMut[t][mut]['cells']
                    diploidPop = 2.0*self.popSize[t]
                    vaf = cells/diploidPop
                    cellsWithMut[t][mut]['rawVAF']=vaf
                except ZeroDivisionError:
                    vaf = 0.0
                    cellsWithMut[t][mut]['rawVAF']=vaf

                if cellsWithMut[t][mut]['rawVAF'] > 1.0:
                    print("DEBUG time: %s"%(t))
                    print("DEBUG Cells: %s"%(cells))
                    print("DEBUG VAF: %s"%(vaf))
                    print("DEBUG MUT: %s"%(mut))
                    print("Population: %s"%(self.popSize[t]))
                    print(self.popsWithMutsAtTimepoint[t][mut])
                    print(cellsWithMut[t][mut])
                assert cellsWithMut[t][mut]['rawVAF'] <= 1.0, "VAF Greater than 1.0 error: %s"%(cellsWithMut[t][mut]['rawVAF'])

        self.VAFs = cellsWithMut
        if(Options.verbose):
            print("Calculating VAFs: %.2f seconds"%(time.clock()-start_time))

        start_time = time.clock()
        self._getCorrectedVAF(Options)
        if(Options.verbose):
            print("Sampling and correcting VAFs: %.2f seconds"%(time.clock()-start_time))

    def _getCorrectedVAF(self, Options):
        '''
        Gets a corrected VAF based on the rawVAF using either a poisson distribution w/
        binomial or a gamma fit dist if options are selected
        :param Options: Command arguments
        :return: None. Sets values in a dictionary.
        '''
        # TODO need to add GAMMA distribution fits from a file
        if Options.depthFile!=None:
            pass
        else:
            for t in self.VAFs:
                for mut in self.VAFs[t]:
                    vaf = self.VAFs[t][mut]['rawVAF']
                    try:
                        self.VAFs[t][mut]['reads'] = binom.rvs(Options.depth, vaf)
                    except ValueError: # TODO vaf shouldn't ever be more than 1.0 but it is.
                        self.VAFs[t][mut]['reads'] = binom.rvs(Options.depth, 1.0)

                    self.VAFs[t][mut]['depth'] = float(poisson.rvs(Options.depth))
                    # print("%s reads, %s depth"%(fi,Di))
                    self.VAFs[t][mut]['correctedVAF'] = self.VAFs[t][mut]['reads']/self.VAFs[t][mut]['depth']

def BuildOutputTable(Options, data, snpeff, refGenome):
    geneLocs = GeneLocs(Options, snpeff)
    pickle.dump(geneLocs, open('file.p', 'wb'))  # TODO Delete when done with dev
    geneLocs = pickle.load(open('file.p', 'rb'))  # TODO Delete when done with dev

    fastaRec = FastaRecords(geneLocs, Options, snpeff, refGenome)
    pickle.dump(fastaRec, open('file2.p', 'wb'))  # TODO Delete when done with dev
    fastaRec = pickle.load(open('file2.p', 'rb'))  # TODO Delete when done with dev
    revDict = {'t': 'a', 'a': 't', 'g': 'c', 'c': 'g'}

    outline = []
    for i, sim in enumerate(data):
        if Options.tmpts==None: # Means we are getting every tmpt
            for k, t in enumerate(sim.tmpts):
                for v, mut in enumerate(sim.mutsAtTimepoint[k]):
                    # Sample    Tmpt    Chrom   Pos Ref Alt Gene    Induction_t Trinuc TrinucTrue  RawVAF  CorrVAF Reads Depth
                    mutInfo = mut.split('.')
                    ref = mutInfo[3].split('[')[1].split('>')[0]
                    alt = mutInfo[3].split('[')[1].split('>')[1].split(']')[0]
                    vafDat = sim.VAFs[k][mut]
                    val = [os.path.basename(sim.data).replace('.csv',''), t, mutInfo[4], mutInfo[5], ref, alt, mutInfo[2], mutInfo[0], mutInfo[3][0] + alt + mutInfo[3][6], mutInfo[3][0] + alt + mutInfo[3][6], vafDat['rawVAF'], vafDat['correctedVAF'], vafDat['reads'], vafDat['depth']]
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

def BuildOutputCloneData(Options, data):
    outline = []
    for i, sim in enumerate(data):
        for clone in sim.parsedData:
            line = [os.path.basename(sim.data).replace('.csv',''), sim.parsedData[clone]['CloneID'], sim.parsedData[clone]['ParentID'], len(sim.parsedData[clone]['Genome'].split(';')), '\t'.join(sim.parsedData[clone]['Pops'])]
            outline.append('\t'.join([str(item) for item in line]))

    with open(Options.output + "mutations.clonal.txt", "w") as outCounts:
        outCounts.write('\t'.join(["Sample","CloneID", "ParentID", "NumMuts",'\t'.join([str(i) for i, v in enumerate('\t'.join(sim.parsedData[0]['Pops']))])]) + "\n")
        outCounts.write("\n".join(outline))

def AnnotateVCF(Options, snpeff, inFile):
    print("Annotating mutations...")
    cmd = ["java", "-jar", snpeff['snpeff'], "ann", "-canon", Options.genome, inFile]
    locInfo = os.popen(' '.join(cmd)).read()
    with open(inFile.replace('.vcf','.ann.vcf'), 'w') as outfile:
        outfile.write(locInfo)
    print("Annotation Complete...")
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

    print("Processing GattacaExample.Gattaca Outputs.")

    if os.path.isdir(os.path.expanduser(Options.output))==False:
        os.mkdir(os.path.expanduser(Options.output))

    data = []
    if os.path.isfile(os.path.expanduser(Options.inFile)):
        data = Data(Options)
    elif os.path.isdir(os.path.expanduser(Options.inFile)):
        allFiles = glob.glob(Options.inFile + "*.csv")
        Options2 = Options
        for i in allFiles:
            Options2.inFile = i
            data.append(Data(Options2))
    else:
        sys.exit("No input file(s) found.")
    #
    #
    # # # TODO build endpoint table or table with all mutations.
    # pickle.dump(data, open('postProcess.p', 'wb')) # TODO Delete when done with dev
    # data = pickle.load(open('postProcess.p', 'rb')) # TODO Delete when done with dev
    BuildOutputCloneData(Options, data)
    BuildOutputTable(Options, data, snpeff, refGenome)

    if Options.plots:
        build = Plots(Options,data)
        # build.mutContextMovie(Options, data)
        # build.mutVAF(Options, data)
        build.oneOverF(Options, data)



