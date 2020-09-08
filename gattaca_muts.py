import sys
import numpy as np

class Mutations:

    def __init__(self, geneData, mutRate):
        self.genes = self._getGenes(geneData)
        self.lawrenceRates = self._getLawrenceRates()
        self.adjRates = None
        self.geneMuts = self._getUserMutRate(geneData, mutRate)

    def _getGenes(self, geneData):
        return([gene for gene in geneData])

    def _getLawrenceRates(self):
        with open("./data/Lawrence_et_al_2013.txt", "r") as inFile:
            vals = inFile.readlines()
        out = dict()
        for line in vals:
            if line.startswith("gene")==False:
                out.update({ line.replace('\n','').split('\t')[0] : float(line.replace('\n','').split('\t')[1]) })
        return(out)

    def _getUserMutRate(self, geneData, mutRate):
        rates = []
        for gene in self.genes:
            try:
                rates.append(self.lawrenceRates[gene])
            except KeyError:
                print("WARNING: Unable to find mutation rate for %s using user defined average."%gene)
                rates.append(mutRate)

        x = np.mean(rates)
        factor = mutRate/x

        userRates = np.asarray(rates)*factor
        self.adjRates = {}
        for i, gene in enumerate(self.genes):
            self.adjRates.update({gene : list(userRates)[i]})

        expectMuts = {}
        for i, gene in enumerate(self.genes):
            expectMuts.update({ gene : (geneData[gene]['end']-geneData[gene]['start'])*userRates[i] })

        return(expectMuts)

