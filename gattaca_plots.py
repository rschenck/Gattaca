import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import sys
from scipy import stats

class MutPropPlots:

    def __init__(self, mutProp, triNucs, triNucNumbers, contextorders, muts, geneLocs):
        self._96mutContexts(mutProp)
        self._32possibleTrinucs(triNucs)
        # self._plotGeneNumbers(triNucNumbers, contextorders)
        # self._plotMutationRates(muts, geneLocs)

    def _96mutContexts(self, mutProp):
        barlist = plt.bar(mutProp.keys(), mutProp.values(), align='center', alpha=0.7)
        plt.xticks(rotation=90)
        self._set_size(10,2)
        plt.tick_params(axis='x', labelsize=8)
        plt.ylim(top=self._getYaxis(mutProp.values()))
        plt.ylabel("Proportion of Mutation")
        plt.xlabel("")
        color = dict()
        theColors = ['blue','black','red','grey','green','pink']
        for i, k in enumerate(theColors):
            color.update({i:k})
        col=0
        for i in range(0,96):
            if i%16==0 and i!=0:
                col+=1
            barlist[i].set_color(color[col])
        # plt.show()
        plt.savefig('imgs/96plot.gattaca.pdf', bbox_inches='tight')
        plt.close()

        # Get one that shows uniform probabilities
        barlist = plt.bar(mutProp.keys(), [1/96. for item in range(0,96)], align='center', alpha=0.7)
        plt.xticks(rotation=90)
        self._set_size(10, 2)
        plt.tick_params(axis='x', labelsize=8)
        plt.ylim(top=self._getYaxis([1/96. for item in range(0,96)]))
        plt.ylabel("Proportion of Mutation")
        plt.xlabel("")
        color = dict()
        theColors = ['blue', 'black', 'red', 'grey', 'green', 'pink']
        for i, k in enumerate(theColors):
            color.update({i: k})
        col = 0
        for i in range(0, 96):
            if i % 16 == 0 and i != 0:
                col += 1
            barlist[i].set_color(color[col])
        # plt.show()
        plt.savefig('imgs/96BlankPlot.gattaca.pdf', bbox_inches='tight')
        plt.close()

        # Single set for a figure: (GCG example)
        # ex = dict()
        # for val in mutProp:
        #     if val == "GCG":
        #         print(val)


    def _32possibleTrinucs(self, triNucs):
        # Get one that shows uniform probabilities
        o = sorted(triNucs.keys())
        x = []
        y = []
        for i in o:
            x.append(i)
            y.append(triNucs[i])

        barlist = plt.bar(x, y, align='center', alpha=0.7)
        # barlist = plt.bar(triNucs.keys(), [1 / 32. for item in range(0, 96)], align='center', alpha=0.7)
        plt.xticks(rotation=90)
        self._set_size(10, 2)
        plt.tick_params(axis='x', labelsize=8)
        plt.ylim(top=self._getYaxis([1 / 96. for item in range(0, 96)]))
        plt.ylabel("Proportion of Mutation")
        plt.xlabel("")
        color = dict()
        theColors = ['blue', 'black', 'red', 'grey', 'green', 'pink', 'purple', 'aqua']
        for i, k in enumerate(theColors):
            color.update({i: k})
        col = 0
        for i in range(0, 32):
            if i % 4 == 0 and i != 0:
                col += 1
            barlist[i].set_color(color[col])
        # plt.show()
        plt.savefig('imgs/32contexts.gattaca.pdf', bbox_inches='tight')
        plt.close()

        barlist = plt.bar(x, [1/32. for i in y], align='center', alpha=0.7)
        # barlist = plt.bar(triNucs.keys(), [1 / 32. for item in range(0, 96)], align='center', alpha=0.7)
        plt.xticks(rotation=90)
        self._set_size(10, 2)
        plt.tick_params(axis='x', labelsize=8)
        plt.ylim(top=self._getYaxis([1 / 96. for item in range(0, 96)]))
        plt.ylabel("Proportion of Mutation")
        plt.xlabel("")
        color = dict()
        theColors = ['blue', 'black', 'red', 'grey', 'green', 'pink', 'purple', 'aqua']
        for i, k in enumerate(theColors):
            color.update({i: k})
        col = 0
        for i in range(0, 32):
            if i % 4 == 0 and i != 0:
                col += 1
            barlist[i].set_color(color[col])
        # plt.show()
        plt.savefig('imgs/32contextsBlank.gattaca.pdf', bbox_inches='tight')
        plt.close()

    def _set_size(self, w,h, ax=None):
        """ w, h: width, height in inches """
        if not ax: ax=plt.gca()
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        figw = float(w)/(r-l)
        figh = float(h)/(t-b)
        ax.figure.set_size_inches(figw, figh)

    def _getYaxis(self, values):
        max = np.max([item for item in values])
        if max > 0.8:
            return(1)
        else:
            return(max+0.15)

    def _plotGeneNumbers(self, g, ordered):
        genes = dict.fromkeys(g)
        for i in g:
            tri = OrderedDict.fromkeys(g[i]) # get a dictionary of the trinucleotide keys
            for item in tri:
                tri[item] = len(g[i][item])/float(sum([len(g[i][k]) for k in g[i]]))
            genes[i] = tri

        genes = g.keys()
        tri = OrderedDict.fromkeys(ordered)  # get a dictionary of the trinucleotide keys
        for i in tri:
            tri[i]=[]

        notPlotted = []
        for gene in genes:
            for item in tri:
                theSum = sum( [len(g[gene]['pos'][k]) for k in g[gene]['pos']] )
                if theSum!=0:
                    tri[item].append( len(g[gene]['pos'][item]) / theSum )
                else:
                    notPlotted.append(gene)
        print(notPlotted)

        colors = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059","#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87","#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80","#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100","#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F","#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09","#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66","#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C","#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81","#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00","#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700","#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329","#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]

        fig, ax = plt.subplots()
        ax.bar(genes, tri[ordered[0]], label=ordered[0], color=colors[0])
        for i, k in enumerate(ordered):
            if i >0:
                sumBeneath = np.zeros(len(genes))
                for v in range(0,len(ordered)):
                    if v<i:
                        sumBeneath+=np.asarray(tri[ordered[v]])
                ax.bar(genes, tri[k],bottom=sumBeneath, color=colors[i], label=ordered[i])

        self._set_size(10, 3)
        ax.set_ylabel("Fraction of Trinucleotides")
        ax.set_xlabel("Gene")
        plt.savefig('imgs/geneTrinucsExample.gattaca.pdf', bbox_inches='tight')
        plt.close()

    def _plotMutationRates(self, muts, genes):
        lengths = []
        lawrenceRates = []
        adjRates = []
        for gene in muts.genes:
            lengths.append((genes.data[gene]['end']-genes.data[gene]['start']))
            lawrenceRates.append(muts.lawrenceRates[gene])
            adjRates.append(muts.adjRates[gene])


        plt.plot(lengths, lawrenceRates, 'go', label="Mutation Rate")
        plt.plot(lengths, adjRates, 'bo', label="Mutation Rate")
        plt.xlabel("Length of gene")
        plt.ylabel("Mutation rate")
        plt.savefig('imgs/mutationRates.gattaca.pdf', bbox_inches='tight')
        plt.close()


