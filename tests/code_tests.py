'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''

from gattaca_classes import GeneLocs

class Options:
    def __init__(self):
        self.geneList = './badGenes.txt'
        self.genome = 'GRCh37.75'

def testSnpEff():
    options = Options()
    geneLocs = GeneLocs(options ,{'snpeff':'/Users/rschenck/Desktop/BioinformaticsTools/snpEff/snpEff.jar'})
    assert len(geneLocs.data)==1, "Improper handling of bad genes."

if __name__=="__main__":
    testSnpEff()