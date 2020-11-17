'''
@author: Ryan Schenck, ryan.schenck@univ.ox.ac.uk
'''
import sys, os
from collections import OrderedDict

class JavaCode:

    def __init__(self, geneList, contextProbs, expectedMuts, triNucsProb, triNucsGenePos, outDir):
        self.orderedVals = self._getProperOrder96(contextProbs, triNucsProb)
        self.lines = self.getJavaCode(geneList, contextProbs, expectedMuts, triNucsProb, triNucsGenePos, outDir)
        self.writeJavaClass(outDir)

    def getJavaCode(self, geneList, contextProbs, expectedMuts, triNucsProb, triNucsGenePos, outDir):
        self._getGenePositions(triNucsGenePos, triNucsProb, geneList, outDir)
        if outDir[len(outDir)-1]!='/':
            outDir=outDir+"/"
        packageName=(outDir.split('/')[len(outDir.split('/'))-2])

        head = "/**\n * Created by Gattaca (Ryan O. Schenck).\n */\n\npackage " + packageName + ";\n\nimport HAL.Tools.FileIO;\nimport cern.jet.random.Poisson;\nimport cern.jet.random.engine.DRand;\nimport cern.jet.random.engine.RandomEngine;\nimport HAL.Tools.PhylogenyTracker.Genome;\nimport HAL.Tools.MultinomialCalc;\nimport HAL.Rand;\n\nimport java.util.ArrayList;\n\n"
        part1 = "public class Gattaca extends Genome<Gattaca> {\n"
        part2 = "\tprivate static final int genomeComponents=%s;"%(len(geneList))
        part3 = "\tprivate static final String[] geneNames=new String[]{%s};"%(','.join(['"'+item+'"' for item in geneList]))
        part3b = "\tprivate static final String[] triNucs=new String[]{%s};"%(','.join(['"'+item+'"' for item in triNucsProb]))
        part3c = "\tprivate static final double[] triNucProbs=new double[]{%s};"%(','.join([repr(triNucsProb[item]) for item in triNucsProb]))
        part3d = "\tprivate static final String[] chrom=new String[]{%s};"%( ','.join(['"'+triNucsGenePos[gene]['chr']+'"' for gene in geneList]) )
        part4 = "\tprivate static final String[] mutTypes=new String[]{%s};"%(','.join([repr(item) for item in self.orderedVals[0]]))
        part4b = "\tprivate static final double[] mutTypeProb=new double[]{%s};"%(','.join([repr(item) for item in self.orderedVals[1]]))
        part4 = part4.replace("'",'"')
        part5 = "\tprivate static double[] expectedMuts=new double[]{%s};"%(','.join([repr(expectedMuts[item]) for item in geneList]))

        vars = "\n\tprivate static final String BaseIndexFile= System.getProperty(\"user.dir\") + \"/src/triNucsPos.csv\";\n\tprivate static final long[][][] basePositions=ParseTriNucFile();\n\tprivate Poisson[] PoissonDists;\n\tprivate RandomEngine RNEngine;\n\tprivate Rand RN;\n\tprivate MultinomialCalc picker;\n\n\tString PrivateGenome;\n\tdouble h;\n\tdouble s;\n\tdouble v;\n\tint[] triNucMuts;\n\tstatic double[] theMutOptions=new double[3];\n\tstatic double[] mutOptionsProbs=new double[3];\n\tint[] mutTypeHit;\n\tint mutIdx;\n\tString mutChosen;\n\tString mutChrom;\n\tString mutGene;\n\tlong mutPos;\n\tint mutations;"

        extraFuncs = "private int NonZeroIdx(int[] arr){\n\t\tint idx=199;\n\t\tint sum=0;\n\t\tfor (int i = 0; i < arr.length; i++) {\n\t\t\tif(arr[i]!=0){\n\t\t\t\tidx=i;\n\t\t\t}\n\t\t\tsum+=arr[i];\n\t\t}\n\t\tif(sum>1 | idx==199){\n\t\t\tSystem.out.println(arr[0] + \"\t\" + arr[1] + \"\t\" + arr[2]);\n\t\t\tSystem.out.println(idx);\n\t\t\tthrow new IllegalStateException(\"Multiple mutation types chosen: \"+sum);\n\t\t}\n\t\treturn(idx);\n\t}\n\n\tprivate void BuildPoissons(){\n\t\tPoissonDists = new Poisson[expectedMuts.length];\n\t\tfor (int i = 0; i < expectedMuts.length; i++) {\n\t\t\tPoisson poisson_dist = new Poisson(expectedMuts[i], RNEngine);\n\t\t\tPoissonDists[i] = poisson_dist;\n\t\t}\n\t}\n\n\tprivate static void NormalizeMutTypeOptions(){\n\t\tdouble sum = sumArray(theMutOptions);\n\t\tfor (int i = 0; i < theMutOptions.length; i++) {\n\t\t\tmutOptionsProbs[i]=theMutOptions[i]/sum;\n\t\t}\n\t\tdouble check=1;\n\t\tfor (int i = 0; i < mutOptionsProbs.length; i++) {\n\t\t\tcheck-=mutOptionsProbs[i];\n\t\t}\n\n\t\tmutOptionsProbs=FixSumToOne(check, mutOptionsProbs);\n\t}\n\n\tprivate static double sumArray(double[] vals){\n\t\tdouble val=0;\n\t\tfor (int i = 0; i < vals.length; i++) {\n\t\t\tval+=vals[i];\n\t\t}\n\t\treturn(val);\n\t}\n\n\tprivate static double[] FixSumToOne(double check, double[] arr){\n\t\tif(check!=0){\n\t\t\tfor (int i = 0; i < arr.length; i++) {\n\t\t\t\tif(check<0){\n\t\t\t\t\tarr[i]-=check/arr.length;\n\t\t\t\t} else if(check>0) {\n\t\t\t\t\tarr[i]+=check/arr.length;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\treturn(arr);\n\t}\n\n\tpublic void GattacaMultinomial(double[] probabilities, int n, int[] ret) {\n\t\tdouble probRemaining = 1;\n\t\tfor (int i = 0; i < probabilities.length; i++) {\n\t\t\tdouble prob=probabilities[i];\n\t\t\tif(n==0){\n\t\t\t\tret[i]=0;\n\t\t\t\treturn;\n\t\t\t}\n\t\t\tif(probRemaining-prob==0){\n\t\t\t\tret[i]=n;\n\t\t\t\treturn;\n\t\t\t}\n\t\t\tint popSelected=picker.Binomial(n,prob/ probRemaining);\n\t\t\tprobRemaining -=prob;\n\t\t\tn-=popSelected;\n\t\t\tret[i]=popSelected;\n\t\t}\n\t}\n\n\t// Parses Base Mutation Information\n\tprivate static long[][][] ParseTriNucFile(){\n\t\tFileIO reader = new FileIO(BaseIndexFile, \"r\");\n\t\tArrayList<long[]> data = new ArrayList<> (reader.ReadLongs(\",\"));\n\t\tlong[][][] BaseIndexes = new long[data.size()/32][32][];\n\t\tfor (int i = 0; i < data.size(); i+=32) {\n\t\t\tBaseIndexes[i/32][0] = data.get(i+0);\n\t\t\tBaseIndexes[i/32][1] = data.get(i+1);\n\t\t\tBaseIndexes[i/32][2] = data.get(i+2);\n\t\t\tBaseIndexes[i/32][3] = data.get(i+3);\n\t\t\tBaseIndexes[i/32][4] = data.get(i+4);\n\t\t\tBaseIndexes[i/32][5] = data.get(i+5);\n\t\t\tBaseIndexes[i/32][6] = data.get(i+6);\n\t\t\tBaseIndexes[i/32][7] = data.get(i+7);\n\t\t\tBaseIndexes[i/32][8] = data.get(i+8);\n\t\t\tBaseIndexes[i/32][9] = data.get(i+9);\n\t\t\tBaseIndexes[i/32][10] = data.get(i+10);\n\t\t\tBaseIndexes[i/32][11] = data.get(i+11);\n\t\t\tBaseIndexes[i/32][12] = data.get(i+12);\n\t\t\tBaseIndexes[i/32][13] = data.get(i+13);\n\t\t\tBaseIndexes[i/32][14] = data.get(i+14);\n\t\t\tBaseIndexes[i/32][15] = data.get(i+15);\n\t\t\tBaseIndexes[i/32][16] = data.get(i+16);\n\t\t\tBaseIndexes[i/32][17] = data.get(i+17);\n\t\t\tBaseIndexes[i/32][18] = data.get(i+18);\n\t\t\tBaseIndexes[i/32][19] = data.get(i+19);\n\t\t\tBaseIndexes[i/32][20] = data.get(i+20);\n\t\t\tBaseIndexes[i/32][21] = data.get(i+21);\n\t\t\tBaseIndexes[i/32][22] = data.get(i+22);\n\t\t\tBaseIndexes[i/32][23] = data.get(i+23);\n\t\t\tBaseIndexes[i/32][24] = data.get(i+24);\n\t\t\tBaseIndexes[i/32][25] = data.get(i+25);\n\t\t\tBaseIndexes[i/32][26] = data.get(i+26);\n\t\t\tBaseIndexes[i/32][27] = data.get(i+27);\n\t\t\tBaseIndexes[i/32][28] = data.get(i+28);\n\t\t\tBaseIndexes[i/32][29] = data.get(i+29);\n\t\t\tBaseIndexes[i/32][30] = data.get(i+30);\n\t\t\tBaseIndexes[i/32][31] = data.get(i+31);\n\t\t}\n\t\treturn BaseIndexes;\n\t}"

        end = "}"
        return("\n".join([head,part1,part2,part3,part5,part3b,part3c,part3d,part4,part4b,vars, self._constructor(), self._mutationFunction(triNucsProb), extraFuncs, end]))

    def writeJavaClass(self, outDir):
        with open(outDir + 'Gattaca.java', 'w') as outfile:
            outfile.write(self.lines)

    def _mutationFunction(self, triNucsProb):
        ret='\tpublic Gattaca _RunPossibleMutation(double time) {\n\t\tStringBuilder MutsObtained = new StringBuilder();\n\t\ttriNucMuts=new int[triNucProbs.length]; // Must reset to zeroes\n\t\tint mutsHappened=0;\n\t\tfor (int j = 0; j < expectedMuts.length; j++) {\n\t\t\tmutations = PoissonDists[j].nextInt(); // How many mutations for this gene.\n\t\t\tif (mutations>0) {\n\t\t\t\tmutsHappened+=mutations;\n\t\t\t\ttriNucMuts = new int[triNucProbs.length];\n\t\t\t\tGattacaMultinomial(triNucProbs, mutations, triNucMuts); // Step One acquire mutation trinucleotide based on 32 mutation types\n\n\t\t\t\tfor (int tri = 0; tri < triNucMuts.length; tri++) {\n\t\t\t\t\tfor (int hit = 0; hit < triNucMuts[tri]; hit++) { // Step Two acquire mutation type based on probability of context transition/transversion to get the 96 mutation types.\n\t\t\t\t\t\t\ttheMutOptions[0] = mutTypeProb[(tri * 3) + 0];\n\t\t\t\t\t\t\ttheMutOptions[1] = mutTypeProb[(tri * 3) + 1];\n\t\t\t\t\t\t\ttheMutOptions[2] = mutTypeProb[(tri * 3) + 2];\n\t\t\t\t\t\t\tNormalizeMutTypeOptions();\n\t\t\t\t\t\t\tmutTypeHit = new int[3];\n\t\t\t\t\t\t\tGattacaMultinomial(mutOptionsProbs, 1, mutTypeHit);\n\t\t\t\t\t\t\tmutIdx = NonZeroIdx(mutTypeHit);\n\t\t\t\t\t\t\tmutChosen = mutTypes[(tri*3)+mutIdx]; // Which mutation occurs in bases.\n\t\t\t\t\t\t\tmutPos = basePositions[j][tri][RN.rn.Int(basePositions[j][tri].length)];\n\t\t\t\t\t\t\tmutChrom = chrom[j];\n\t\t\t\t\t\t\tmutGene = geneNames[j];\n\t\t\t\t\t\t\tString mutOut = time + "." + mutGene + "." + mutChosen + "." + mutChrom + "." + mutPos + ";";\n\t\t\t\t\t\t\tMutsObtained.append(mutOut);\n\t\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\n\t\t}\n\n\t\tif(mutsHappened>0){\n\t\t\tString NewPrivateGenome=MutsObtained.toString();\n\t\t\treturn(new Gattaca(this, NewPrivateGenome, RN.Double(),RN.Double(),1,RN));\n\t\t} else {\n\t\t\treturn(this);\n\t\t}\n\t}'

        return(ret)

    def _constructor(self):
        constructor =  "\tpublic Gattaca(Gattaca parent, String PrivateGenome, double h, double s, double v, Rand RN){\n\t\tsuper(parent, false);\n\t\tthis.PrivateGenome = PrivateGenome;\n\t\tthis.RNEngine=new DRand();\n\t\tthis.RN=RN;\n\t\tthis.picker=new MultinomialCalc(RN);\n\t\tthis.BuildPoissons();\n\t\tthis.h=h;\n\t\tthis.s=s;\n\t\tthis.v=v;\n\t}"
        return(constructor)

    def _getProperOrder96(self, contextProbs, triNucsProb):
        ordered = []
        orderedValue = []
        for tri in triNucsProb:

            for triMut in contextProbs:

                if tri == triMut[0].lower()+triMut[2].lower()+triMut[6].lower():
                    ordered.append(triMut)
                    orderedValue.append(contextProbs[triMut])

        return((ordered,orderedValue))

    def _getGenePositions(self, positions, triNucsProb, geneList, outDir):
        '''
        Access positions based on geneList, order in triNucsProb order. This will take the shape of:
        long[][][] where sizes are [gene][32contexts][allPositions] so { { { } }    }
        '''
        # bases = "private static final long[][][] basePositions=new long[][][]{"
        #
        # for j, gene in enumerate(geneList):
        #     outTris = OrderedDict()
        #     for tri in triNucsProb:
        #         thisTri = "{"
        #         for i, pos in enumerate(positions[gene]['pos'][tri]):
        #             if i != len(positions[gene]['pos'][tri])-1:
        #                 thisTri+=repr(pos)+"l,"
        #             else:
        #                 thisTri += repr(pos) + "l}"
        #
        #         outTris.update({tri:thisTri})
        #
        #     geneOut = "\n{\n"
        #     for i, outTri in enumerate(outTris):
        #         if i != len(outTris) - 1:
        #             geneOut += outTris[outTri]+",\n"
        #         else:
        #             geneOut += outTris[outTri]
        #
        #     if j != len(geneList)-1:
        #         geneOut += "\n},"
        #     else:
        #         geneOut += "\n}"
        #
        #     bases += geneOut
        #
        # bases+="\n}"

        with open(outDir + "triNucsPos.csv", "w") as outFile:
            for j, gene in enumerate(geneList):
                print("Writing gene position file: %s" %(gene))

                for tri in triNucsProb:
                    thisTri = []
                    for i, pos in enumerate(positions[gene]['pos'][tri]):
                        thisTri.append(repr(pos))

                    if thisTri == []:
                        print("ERROR No Positions found: %s" % (gene))

                    outFile.write(','.join(thisTri) + '\n')

        # return(bases)

