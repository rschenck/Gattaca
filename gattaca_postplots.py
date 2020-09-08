import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import sys
import scipy.stats as st
import os
import subprocess
import glob
import time
import math

class Plots:

    def __init__(self, Options, data):
        pass

    def setVAFs(self, Options, data):
        if type(data)==type([]):
            for t, tmpt in enumerate(data[0].tmpts): # Looping over timepoints
                self.rawVAFs = []
                self.correctedVAFs = []
                self.depth = []
                for i, val in enumerate(data): # Looping over simulations
                    for k, mut in enumerate(data[i].VAFs[tmpt]): # Looping over mutations
                        # print(data[i].VAFs[tmpt][mut])
                        self.rawVAFs.append(data[i].VAFs[tmpt][mut]['rawVAF'])
                        self.correctedVAFs.append(data[i].VAFs[tmpt][mut]['correctedVAF'])
                        self.depth.append(data[i].VAFs[tmpt][mut]['depth'])

    def mutContextMovie(self, Options, data):
        for t, tmpt in enumerate(data[0].tmpts):
            allMuts = 0
            out = OrderedDict.fromkeys(data[0].contexts)
            for k in out:
                out[k] = 0
            for i, dataClass in enumerate(data):
                for tri in out:
                    allMuts+=dataClass.mutContexts[t][tri]
                    out[tri]+=dataClass.mutContexts[t][tri]
            for k in out:
                try:
                    out[k]=out[k]/allMuts
                except ZeroDivisionError:
                    pass

            plt.figure(figsize=(12,4) , dpi=72)
            barlist = plt.bar(out.keys(), out.values(), align='center', alpha=0.7)
            plt.xticks(rotation=90)
            plt.tick_params(axis='x', labelsize=8)
            plt.ylim(bottom=0, top=getYaxis(out.values()))
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
            plt.savefig(Options.output + '96plot.'+str.zfill(str(t), len(str(len(data[0].tmpts))))+'.png', bbox_inches='tight', dpi=300)
            plt.close()

        pad = len(str(len(data[0].tmpts)))
        # generate_video(Options, pad)
    #     ffmpeg -framerate 8 -i 96plot.%02d.png -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./96mutContextPlot.mp4 -vcodec libx264

    def mutVAF(self, Options, data):
        print("Generating Muation VAF plots.")
        if type(data)==type([]):
            for t, tmpt in enumerate(data[0].tmpts): # Looping over timepoints
                self.rawVAFs = []
                self.correctedVAFs = []
                self.depth = []
                for i, val in enumerate(data):  # Looping over simulations
                    for k, mut in enumerate(data[i].VAFs[tmpt]):  # Looping over mutations
                        # print(data[i].VAFs[tmpt][mut])
                        self.rawVAFs.append(data[i].VAFs[tmpt][mut]['rawVAF'])
                        self.correctedVAFs.append(data[i].VAFs[tmpt][mut]['correctedVAF'])
                        self.depth.append(data[i].VAFs[tmpt][mut]['depth'])

                fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=False, sharex=False)

                ax1.hist(self.rawVAFs, density=True, bins=60, label="Data", color="darkgreen")
                mn, mx = plt.xlim()
                try:
                    kde_xs = np.linspace(mn, max(self.rawVAFs)+0.1, 301)
                except ValueError:
                    kde_xs = np.linspace(mn, mx, 301)
                try:
                    kde = st.gaussian_kde(self.rawVAFs)
                    ax1.plot(kde_xs, kde.pdf(kde_xs), label="PDF", color="black")
                except ValueError:
                    pass
                ax1.legend(loc="upper right")
                ax1.set_title('Raw (True) VAF')
                ax1.set_xlabel("VAF")
                ax1.set_ylabel("Count")

                ax2.hist(self.correctedVAFs, density=True, bins=60, label="Data", color="darkgreen")
                mn, mx = plt.xlim()
                try:
                    kde_xs = np.linspace(mn, max(self.correctedVAFs)+0.1, 301)
                except ValueError:
                    kde_xs = np.linspace(mn, mx, 301)
                try:
                    kde = st.gaussian_kde(self.correctedVAFs)
                    ax2.plot(kde_xs, kde.pdf(kde_xs), label="PDF", color="black")
                except ValueError:
                    pass
                ax2.legend(loc="upper right")
                ax2.set_title('Corrected (Sampled) VAF')
                ax2.set_xlabel("VAF")
                ax2.set_ylabel("Count")

                ax3.hist([val for val in self.correctedVAFs if val>=Options.cutoff], density=True, bins=60, label="Data", color="darkgreen")
                mn, mx = plt.xlim()
                try:
                    kde_xs = np.linspace(mn, max([val for val in self.correctedVAFs if val>=Options.cutoff])+0.1, 301)
                except ValueError:
                    kde_xs = np.linspace(mn, mx, 301)
                try:
                    kde = st.gaussian_kde([val for val in self.correctedVAFs if val>=Options.cutoff])
                    ax3.plot(kde_xs, kde.pdf(kde_xs), label="PDF", color="black")
                except ValueError:
                    pass
                ax3.legend(loc="upper right")
                ax3.set_title('Filtered VAF')
                ax3.set_xlabel("VAF")
                ax3.set_ylabel("Count")

                y = [self.depth[i] for i,val in enumerate(self.rawVAFs) if val<Options.cutoff]
                x = [val for i,val in enumerate(self.rawVAFs) if val<Options.cutoff]
                x2 = [val for i,val in enumerate(self.correctedVAFs) if val>=Options.cutoff and val!=0]
                y2 = [self.depth[i] for i,val in enumerate(self.correctedVAFs) if val>=Options.cutoff and val!=0]

                ax4.scatter(x , y, alpha=0.2, s=1.5, color="darkgreen")
                ax4.scatter(x2, y2, alpha=1.0, s=1.5, color="darkgreen")
                ax4.set_title('Reads for correction')
                ax4.set_xlabel("Variant Reads")
                ax4.set_ylabel("Total Reads")

                set_size(10, 3)
                fig.tight_layout(pad=2.0)

                plt.savefig(Options.output + 'VAFplot.' + str.zfill(str(t), len(str(len(data[0].tmpts)))) + '.png',
                            bbox_inches='tight', dpi=120)
                plt.close()

            pad = len(str(len(data[0].tmpts)))
            # generate_video(Options, pad, fileNames="VAFplot.%0", outputName='VAFplot.mp4')
    #         ffmpeg -framerate 8 -i VAFplot.%02d.png -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./VAFplot.mp4 -vcodec libx264

    def mutationBoxPlots(self):
        pass

    def oneOverF(self, Options, data):
        print("Generating Inverse Frequency plots.")
        if type(data)==type([]):
            t = len(data[0].tmpts) - 1
            tmpt = data[0].tmpts[t]
            self.rawVAFs = []
            self.correctedVAFs = []
            self.depth = []
            for i, val in enumerate(data):  # Looping over simulations
                for k, mut in enumerate(data[i].VAFs[t]):  # Looping over mutations
                    # print(data[i].VAFs[tmpt][mut])
                    self.rawVAFs.append(data[i].VAFs[t][mut]['rawVAF'])
                    self.correctedVAFs.append(data[i].VAFs[t][mut]['correctedVAF'])
                    self.depth.append(data[i].VAFs[t][mut]['depth'])


            # roundedRaw = [np.round(v,3) for v in self.rawVAFs]
            # roundedCorrect = [np.round(v,3) for v in self.correctedVAFs]

            # Ensure cutoffs have hard filtered data already
            raw = [val for val in self.rawVAFs if val >= 0.001]
            corrected = [val for val in self.correctedVAFs if val >= Options.cutoff]

            # Set up the bins
            binsRaw = np.linspace(start=max(self.rawVAFs), stop=0.001, num=1000)
            digitizedRaw = np.digitize(raw, binsRaw)

            binsCorrect = np.linspace(start=max(raw), stop=Options.cutoff, num=1000)
            digitizedCorrect = np.digitize(corrected, binsCorrect)

            max_y = max(max(get_cumsum(binsRaw, digitizedRaw) / 100.), max(get_cumsum(binsCorrect, digitizedCorrect) / 100.))
            max_x = max(max(1. / binsRaw), max(1. / binsCorrect)) # will always be the maximum.
            min_y = 0.
            min_x = min(min(1. / binsRaw), min(1. / binsCorrect))  # will always be the maximum.

            for t, tmpt in enumerate(data[0].tmpts):  # Looping over timepoints
                self.rawVAFs = []
                self.correctedVAFs = []
                self.depth = []
                for i, val in enumerate(data):  # Looping over simulations
                    for k, mut in enumerate(data[i].VAFs[tmpt]):  # Looping over mutations
                        # print(data[i].VAFs[tmpt][mut])
                        self.rawVAFs.append(data[i].VAFs[tmpt][mut]['rawVAF'])
                        self.correctedVAFs.append(data[i].VAFs[tmpt][mut]['correctedVAF'])
                        self.depth.append(data[i].VAFs[tmpt][mut]['depth'])

                # roundedRaw = [np.round(v,3) for v in self.rawVAFs]
                # roundedCorrect = [np.round(v,3) for v in self.correctedVAFs]

                # Ensure cutoffs have hard filtered data already
                raw = [val for val in self.rawVAFs if val>=0.001]
                corrected = [val for val in self.correctedVAFs if val>=Options.cutoff]

                # Set up the bins
                try:
                    binsRaw = np.linspace(start=max(self.rawVAFs), stop=0.001, num=1000)
                    digitizedRaw = np.digitize(raw, binsRaw)

                    binsCorrect = np.linspace(start=max(raw), stop=Options.cutoff, num=1000)
                    digitizedCorrect = np.digitize(corrected, binsCorrect)

                    cumSumRaw=get_cumsum(binsRaw, digitizedRaw)/100.
                    cumSumCorrect=get_cumsum(binsCorrect, digitizedCorrect)/100.

                    binsRawPlot = 1./binsRaw
                    binsCorrectPlot = 1./binsCorrect

                    ticks=[min(binsRawPlot),max(binsRawPlot)]
                    labels=["1/%.4f"%max(binsRaw),"1/%.4f"%min(binsRaw)]
                except ValueError:
                    binsRawPlot = []
                    cumSumRaw = []
                    binsCorrectPlot = []
                    cumSumCorrect = []
                    ticks = []
                    labels = []


                fig, (ax1) = plt.subplots(1, 1)
                ax1.scatter(binsRawPlot, cumSumRaw, marker='x', alpha=1.0, s=2, color="darkblue", label="Uncorrected")
                ax1.scatter(binsCorrectPlot, cumSumCorrect, marker='P', alpha=1.0, s=2, color="darkgreen", label="Corrected")
                ax1.set_title('Inverse VAF')
                ax1.set_xlabel("Inverse Variant Frequency (1/f)")
                ax1.set_ylabel("Cumulative Count (x 10^2)")
                ax1.legend(loc="upper left")
                plt.xticks(ticks=ticks, labels=labels)
                plt.axvline(x=1./0.005, linewidth=1, color="black", linestyle="dashed")

                plt.ylim(bottom=0, top=max_y)
                plt.xlim(left=min_x, right=max_x)
                set_size(4,4)
                fig.tight_layout(pad=2.0)

                plt.savefig(Options.output + 'InverseFreq.' + str.zfill(str(t), len(str(len(data[0].tmpts)))) + '.png',
                            bbox_inches='tight', dpi=120)
                plt.close()

                if t==len(data[0].tmpts) - 1:
                    for i, val in enumerate(cumSumCorrect):
                        print("%s\t%s\tCorrected\tFull"%(binsCorrectPlot[i],cumSumCorrect[i]))
                    for i, val in enumerate(cumSumRaw):
                        print("%s\t%s\tRaw\tFull"%(binsRawPlot[i], cumSumRaw[i]))

        pad = len(str(len(data[0].tmpts)))
        # generate_video(Options, pad, fileNames="InverseFreq.%0", outputName="InverseFreqPlot.mp4")
        # ffmpeg -framerate 8 -i InverseFreq.%02d.png -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./InverseFreqPlot.mp4 -vcodec libx264

def get_cumsum(bins, digitized):
    cumSum = np.zeros(len(bins))
    countRaw = 0
    for i, val in enumerate(bins):
        counted = len(np.where(digitized == i)[0])
        countRaw += counted
        cumSum[i] = countRaw
    return(cumSum)

def generate_video(Options, pad, fileNames="96plot.%0", outputName='96mutContextPlot.mp4'):
    os.chdir(Options.output)
    subprocess.call([
        'ffmpeg', '-framerate', '8', '-i', fileNames + str(pad) + 'd.png', '-r', '30', '-pix_fmt', 'yuv420p',
        '-vf', "\"pad=ceil(iw/2)*2:ceil(ih/2)*2\"", outputName
    ])
    for file_name in glob.glob(Options.output+"96plot.*.png"):
        os.remove(file_name)

def getYaxis(values):
    max = np.max([item for item in values])
    if max > 0.8:
        return(1)
    else:
        return(max+0.15)

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = round_up_to_even(float(w)/(r-l))
    figh = round_up_to_even(float(h)/(t-b))
    ax.figure.set_size_inches(figw, figh)

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2
