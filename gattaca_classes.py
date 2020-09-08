'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''
from collections import OrderedDict
import os, sys
import struct
import gzip
from Bio import SeqIO

class GeneLocs:
    '''
    Gets the bedfile of genes of interest
    '''

    def __init__(self, Options, snpeff):
        self.data = self._getBedFileOfGenes(Options, snpeff)

    def _getBedDict(self, results):
        '''
        Builds a dictionary from snpEff genes2bed

        :param results: results from getBedFileOfGenes
        :return: An ordered dictionary of information
        '''
        geneInfo = OrderedDict()
        d = results.split('\n')[1:] # removes header
        for item in d:
            item = [var for var in item.split('\t')]
            if item!=['']:
                info = item[3].split(';')
                v = {'chr':item[0],'start':int(item[1]),'end':int(item[2]),'strand':info[2],'ensembl':info[1]}
                geneInfo.update({ info[0] : v })
        return(geneInfo)


    def _getBedFileOfGenes(self, Options, snpeff):
        '''
        Constructs a bed file from the list of genes.

        :param Options: Parsed options from command line
        :param snpeff: snpeff path from parsed usr_paths.ini
        :return: None
        '''
        with open(Options.geneList, 'r') as infile:
            genes = infile.readlines()
            genes = [v.replace('\n','') for v in genes]

        print("Getting genomic positions for genes...")

        cmd = ["java","-jar", snpeff['snpeff'], "genes2bed", "-t", Options.genome, ' '.join(genes)]
        # print(' '.join(cmd))
        locInfo = os.popen(' '.join(cmd)).read()

        geneInfo = self._getBedDict(locInfo)

        if len(geneInfo)==0:
            sys.exit("No gene information obtained. Ensure that snpEff path and gene names are correct.")

        notFound = 'Genes not found: '
        i=0
        for gene in genes:
            try:
                t = geneInfo[gene]
            except KeyError as k:
                notFound+=gene + ' '
                i+=1
        if i>0:
            print(notFound)
        else:
            print("Records for gene locations found.")

        return(geneInfo)

class FastaRecords:

    def __init__(self, geneLocs, Options, snpeff, ref):
        self.genomeDir = '/'.join(snpeff['snpeff'].split('/')[:-1 or None]) + "/data/" + Options.genome + "/"
        assert os.path.exists(self.genomeDir), "Directory with reference genome is not found."
        self.geneSeq = self._getFastaFileFromGeneLoc(geneLocs, ref)

    def _getFastaFileFromGeneLoc(self, geneLocs, ref):
        chromsToIter = self._setupGeneListByChromosome(geneLocs)
        print("Extracting Sequences...")
        ret = []
        for chrom in chromsToIter:
            print(chromsToIter[chrom])
            ranges, seqs, fastaid, fullseq = self._readChromsomeFile(chrom, ref)

            for gene in chromsToIter[chrom]:
                geneInfo = geneLocs.data[gene]
                print("Attempting to extract sequence for %s %s:%s-%s"%(gene,geneInfo['chr'],geneInfo['start'],geneInfo['end']))
                notFound = True
                for i,v in enumerate(ranges):
                    if ranges[i][0]<=geneInfo['start'] and ranges[i][1]>=geneInfo['end']: # gene is entirely within one scaffolding
                        geneSeq = self._extractSequence(geneInfo['start'], geneInfo['end'],  ranges[i], seqs[i], fastaid, fullseq, gene)
                        notFound = False
                        ret.append(geneSeq)
                    elif ranges[i][0]<=geneInfo['start'] and abs(geneInfo['end']-ranges[i][1])<=1 and abs(geneInfo['end']-ranges[i][1])>=0:
                        geneSeq = self._extractSequence(geneInfo['start'], geneInfo['end'], ranges[i], seqs[i], fastaid,
                                                        fullseq, gene)
                        notFound = False
                        ret.append(geneSeq)
                if notFound:
                    print("WARNING: Gene sequence not recovered from sequences (%s)"%(gene))
        return(ret)

    def _extractSequence(self, start, end,  ranges, seqs, fastaID, fullSeq, gene,  extend=1):
        '''
        Extracts the sequence from the sequence record

        :param start: idx of sequence
        :param stop: idx of end of sequence
        :param range: tuple of sequence coords
        :param seqs: raw sequence of chromosome
        :param extend: extend bases by this many (1 for trinucleotide context)
        :return: values
        '''
        trueStart = start
        trueEnd = end
        geneSeq = seqs[(trueStart-ranges[0]):(ranges[1]-trueEnd)]
        return((geneSeq,trueStart,trueEnd,gene))

    def _readChromsomeFile(self, chr, ref):
        '''
        Retrieves and organizes the chromosome sequence file.

        :param chr: chromosome of interest
        :return: a properly formatted chromosome scaffold
        '''
        chromFile = self.genomeDir + 'sequence.'+chr+'.bin'

        values = []  # Create a list of tuples for the ranges in the chromosome
        seqs = []
        seqname = ""
        fullChrom = ""

        print("Extracting DNA sequence for chromosome %s"%(chr))

        if os.path.exists(chromFile)==False:
            print("Chromosome file not found for %s" % (chr))
            return (values, seqs, seqname, fullChrom)

        chromInfo = []
        with open(chromFile, 'rb') as fd:
            gzip_fd = gzip.GzipFile(fileobj=fd)
            for line in gzip_fd:
                line = line.decode().replace('\n','')
                line=line.split('\t')
                chromInfo.append(line)

        chromInfo = chromInfo[3:]
        values = [] # Create a list of tuples for the ranges in the chromosome
        seqs = []
        for s in chromInfo:
            start = s[3]
            end = s[4]
            seq = s[7]
            values.append((int(start),int(end)))
            seqs.append(seq)

        # Will have to get genome files to complete sequence from assembly.
        fullChrom = ""
        seqname = ""
        # gatherSequence = False
        # with open(ref["ref"], 'rb') as fd:
        #     gzip_fd = gzip.GzipFile(fileobj=fd)
        #     for line in gzip_fd:
        #         line = line.decode()
        #
        #         if (line.startswith(">") and line.split(' ',1)[0].replace('>','')==chr) or gatherSequence==True:
        #             if seqname=="":
        #                 seqname=line.replace('\n','')
        #             gatherSequence = True
        #             if line.startswith(">") and line.split(' ',1)[0].replace('>','')!=chr:
        #                 break
        #             elif seqname!="":
        #                 fullChrom+=line.replace('\n','')

        return(values, seqs, seqname, fullChrom)

    def _setupGeneListByChromosome(self, geneLocs):
        '''
        Organizes genes by which chromosome they are on
        :return: Dictionary chrom:[genes]
        '''
        chroms = list(set([geneLocs.data[gene]['chr'] for gene in geneLocs.data]))
        out = dict.fromkeys(chroms)
        for k in out:
            out[k]=[]
        for chrom in chroms:
            for gene in geneLocs.data:
                if geneLocs.data[gene]['chr']==chrom:
                    out[chrom].append(gene)
        return(out)

class Contexts:

    def __init__(self, Options, geneLocs, fastaRec):
        # self.contextOrder = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]
        self.contextOrder = self._getOrder()
        self.context32Order = ['aca', 'acc', 'acg', 'act', 'ata', 'atc', 'atg', 'att', 'cca', 'ccc', 'ccg', 'cct', 'cta', 'ctc', 'ctg', 'ctt', 'gca', 'gcc', 'gcg', 'gct', 'gta', 'gtc', 'gtg', 'gtt', 'tca', 'tcc', 'tcg', 'tct', 'tta', 'ttc', 'ttg', 'ttt']
        self._getOrder()
        if Options.contextFile=="":
            # self.context=[1/96. for i in range(0,len(self.contextOrder))] # possible contexts, uniform probability
            self.context = self._lawrenceContexts()
        else:
            self.context=self._parseContextFile(Options.contextFile)

        '''Determine the trinucleotide position probability first'''
        self.triNucsProb = self._digestContexts() # should be 32 total. Three bases per

        self.triNucsGenePos = self._splitSeqIntoContexts(fastaRec, geneLocs)

    def _getOrder(self):
        bases = ['A','C','G','T']
        muts = ['[C>A]','[C>G]','[C>T]','[T>A]','[T>C]','[T>G]']
        vals = []
        for pos2 in range(0, 6):
            for pos1 in range(0,4):
                for pos3 in range(0,4):
                        vals.append(bases[pos1]+muts[pos2]+bases[pos3])
        return(vals)

    def _parseContextFile(self, inFile):
        assert os.path.isfile(inFile), "Mutation trinucleotide context file is not found."
        with open(inFile, 'r') as cF:
            lines = cF.readlines()

        lines = OrderedDict([item.replace("\n",'').split('\t') for item in lines])
        linesOut = OrderedDict.fromkeys(self.contextOrder)
        for item in linesOut:
            linesOut[item]=0
        for item in lines:
            try:
                linesOut[item] = float(lines[item])
            except ValueError:
                sys.exit("Unable to convert string to float for mutation contexts.")

        return(linesOut)

    def _splitSeqIntoContexts(self, fastaRec, geneLocs):
        contexts = {}
        for item in fastaRec.geneSeq:
            gene = item[3]
            seq = item[0]
            start = item[1]
            end = item[2]
            chrom = geneLocs.data[gene]['chr']

            values = self._contextDict(seq, start, end)

            contexts.update({gene:{'chr':chrom, 'pos':values}})
        return(contexts)

    def _contextDict(self, seq, start, end):
        '''
        Gets a list of lists where each entry corresponds to a context list of base pair positions of that context.
        :param seq: sequence of gene
        :param start: start of sequence
        :param end: end of the sequence
        :return: Dictionary of context locations (middle base is the base of that mutation type would be mutated)
        '''
        revcomp = self._revComp(seq)

        ret = OrderedDict.fromkeys(self.context32Order)
        for item in ret:
            ret[item] = []
        for i in range(1,len(seq)-1):
            if seq[i] in ['t','c']:
                trinuc = seq[i-1]+seq[i]+seq[i+1]
            else:
                trinuc = self._revComp(seq[i-1]+seq[i]+seq[i+1])
            ret[trinuc].append(i+(start+1))
        return(ret)

    def _revComp(self, seq):
        revDict = {'t':'a','a':'t','g':'c','c':'g'}
        seqRev = ''.join([revDict[item] for item in seq])
        return(seqRev)

    def _digestContexts(self):
        '''
        For the 32 possible base trinucleotides. This calculates their probability. This is before determining which mutation it becomes.
        :return: 32 probabilities for these trinucleotide choices.
        '''
        values32Ordered = OrderedDict.fromkeys(self.context32Order)
        for entry in values32Ordered:
            values32Ordered[entry]=0
        for c in self.contextOrder:
            match = c[0].lower() + c.split('[')[1].split('>')[0].lower() + c[6].lower()
            values32Ordered[match] += self.context[c]

        assert abs(1-sum(values32Ordered.values())) < 0.01, "Trinucleotide contexts do not sum to one"

        epsilon = 1-sum(values32Ordered.values())

        values32Ordered[entry]+=epsilon
        # print(sum(values32Ordered.values())) # TODO this should always be 1.0 uncomment to check

        return(values32Ordered)

    def _lawrenceContexts(self):
        with open("data/Lawrence_et_al_2013_S3.txt", "r") as inFile:
            lines = [line.replace('\n','').split('\t') for line in inFile.readlines()]

        header=lines[0]
        data=lines[1:]

        keys = []
        for col in header:
            assert col[0]==col[5], "Unknown bases"
            assert col[2]==col[7], "Unknown bases"
            keys.append(col[0]+"["+col[1] + ">" + col[6] + "]" +col[2])

        ready = []
        for sam in data:
            ready.append({key:int(value) for key, value in zip(keys, sam)})

        probs = dict.fromkeys(self.contextOrder)
        for prob in probs:
            probs[prob]=0

        for sample in ready:
            for key in sample:
                try:
                    probs[key]+=sample[key]
                except KeyError:
                    seq = key[0] + key[2] + key[4] + key[6]
                    seqRev = self._revComp(seq.lower()).upper()
                    probs[seqRev[0]+"["+seqRev[1] + ">" + seqRev[2] + "]" +seqRev[3]] += sample[key]

        d = sum(probs.values())
        for key in probs:
            probs[key]=probs[key]/float(d)

        return(probs)
