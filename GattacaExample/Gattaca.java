/**
 * Created by Gattaca (Ryan O. Schenck).
 */

package GattacaExample;

import HAL.Tools.FileIO;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;
import HAL.Tools.PhylogenyTracker.Genome;
import HAL.Tools.MultinomialCalc;
import HAL.Rand;

import java.util.ArrayList;


public class Gattaca extends Genome<Gattaca> {

	private static final int genomeComponents=3;
	private static final String[] geneNames=new String[]{"PIK3CA","KRAS","TP53"};
	private static double[] expectedMuts=new double[]{0.00036969310843373496,0.0001583488,5.536848192771085e-05};
	private static final String[] triNucs=new String[]{"aca","acc","acg","act","ata","atc","atg","att","cca","ccc","ccg","cct","cta","ctc","ctg","ctt","gca","gcc","gcg","gct","gta","gtc","gtg","gtt","tca","tcc","tcg","tct","tta","ttc","ttg","ttt"};
	private static final double[] triNucProbs=new double[]{0.0209303274009191,0.033557599861572876,0.107302079398366,0.02651118778570226,0.038192694284793736,0.00567481172225612,0.013666986802849355,0.01912594074717787,0.02110999360892538,0.01782595926053107,0.06105193293527706,0.017014238588760208,0.0077541000759324494,0.011900540850508678,0.01380158840616718,0.01567323350041558,0.01353409731499741,0.04834403973688821,0.10864636613895277,0.02029160649891354,0.03453495875709209,0.005725065352770743,0.02436827964140632,0.027756541352567873,0.04379543141661486,0.05025955142286499,0.11449695582541337,0.04084707012907112,0.015020886337026971,0.0050083076272755,0.003957482948775858,0.012320144269213486};
	private static final String[] chrom=new String[]{"3","12","17"};
	private static final String[] mutTypes=new String[]{"A[C>A]A","A[C>G]A","A[C>T]A","A[C>A]C","A[C>G]C","A[C>T]C","A[C>A]G","A[C>G]G","A[C>T]G","A[C>A]T","A[C>G]T","A[C>T]T","A[T>A]A","A[T>C]A","A[T>G]A","A[T>A]C","A[T>C]C","A[T>G]C","A[T>A]G","A[T>C]G","A[T>G]G","A[T>A]T","A[T>C]T","A[T>G]T","C[C>A]A","C[C>G]A","C[C>T]A","C[C>A]C","C[C>G]C","C[C>T]C","C[C>A]G","C[C>G]G","C[C>T]G","C[C>A]T","C[C>G]T","C[C>T]T","C[T>A]A","C[T>C]A","C[T>G]A","C[T>A]C","C[T>C]C","C[T>G]C","C[T>A]G","C[T>C]G","C[T>G]G","C[T>A]T","C[T>C]T","C[T>G]T","G[C>A]A","G[C>G]A","G[C>T]A","G[C>A]C","G[C>G]C","G[C>T]C","G[C>A]G","G[C>G]G","G[C>T]G","G[C>A]T","G[C>G]T","G[C>T]T","G[T>A]A","G[T>C]A","G[T>G]A","G[T>A]C","G[T>C]C","G[T>G]C","G[T>A]G","G[T>C]G","G[T>G]G","G[T>A]T","G[T>C]T","G[T>G]T","T[C>A]A","T[C>G]A","T[C>T]A","T[C>A]C","T[C>G]C","T[C>T]C","T[C>A]G","T[C>G]G","T[C>T]G","T[C>A]T","T[C>G]T","T[C>T]T","T[T>A]A","T[T>C]A","T[T>G]A","T[T>A]C","T[T>C]C","T[T>G]C","T[T>A]G","T[T>C]G","T[T>G]G","T[T>A]T","T[T>C]T","T[T>G]T"};
	private static final double[] mutTypeProb=new double[]{0.00488140223956781,0.000804488969084387,0.0152444361922669,0.00802416561805448,0.0181805665776302,0.0073528676658882,0.0133644382724535,0.0114450107078053,0.0824926304181072,0.00956939639731738,0.00474780104243128,0.0121939903459536,0.0102587510850119,0.0252160659511617,0.00271787724862014,0.00110455083837638,0.00457026088387974,0.0,0.0103016667830403,0.00252918215757991,0.000836137862229144,0.00968707322373083,0.00790256735010489,0.00153630017334215,0.010403970413985,0.0018959383088408,0.00881008488609958,0.0123886197963046,0.00205190449619982,0.00338543496802665,0.00308683612332815,0.00377738408789601,0.0541877127240529,0.00183215089781896,0.00244991081553465,0.0127321768754066,0.0,0.00411791571968485,0.0036361843562476,0.00034046654054141,0.00661217252976477,0.0049479017802025,0.00179953602192932,0.00656216208423937,0.00543989029999849,0.00213623486690682,0.00827274854035444,0.00526425009315432,0.00568534628997891,0.00192009324350455,0.00592865778151395,0.00283731160759142,0.00224094228546769,0.0432657858438291,0.010598098950561,0.00182259862606548,0.0962256685623263,0.00109877643596122,0.00679720084580522,0.0123956292171471,0.0145200481728849,0.0107594463752022,0.00925546420900499,0.000560791143230503,0.00516427420954024,0.0,0.00770518639753221,0.0122238575519231,0.00443923569195101,0.00114796903591997,0.0199230273354852,0.0066855449811627,0.0113538051643074,0.00959209784421696,0.0228495284080905,0.00468555342226589,0.02632097722303,0.0192530207775691,0.0122617523579577,0.00129719766875367,0.100938005798702,0.00274229991293782,0.0170005081133654,0.0211042621027679,0.00400406965851523,0.0073075867557355,0.00370922992277624,0.0,0.00288669036059096,0.00212161726668454,0.00130389531129614,0.000441769793851428,0.00221181784362829,0.000923547914709982,0.00439671903461965,0.00699987731988341};

	private static final String BaseIndexFile= System.getProperty("user.dir") + "/src/triNucsPos.csv";
	private static final long[][][] basePositions=ParseTriNucFile();
	private Poisson[] PoissonDists;
	private RandomEngine RNEngine;
	private Rand RN;
	private MultinomialCalc picker;

	String PrivateGenome;
	double h;
	double s;
	double v;
	int[] triNucMuts;
	static double[] theMutOptions=new double[3];
	static double[] mutOptionsProbs=new double[3];
	int[] mutTypeHit;
	int mutIdx;
	String mutChosen;
	String mutChrom;
	String mutGene;
	long mutPos;
	int mutations;
	public Gattaca(Gattaca parent, String PrivateGenome, double h, double s, double v, Rand RN){
		super(parent, false);
		this.PrivateGenome = PrivateGenome;
		this.RNEngine=new DRand();
		this.RN=RN;
		this.picker=new MultinomialCalc(RN);
		this.BuildPoissons();
		this.h=h;
		this.s=s;
		this.v=v;
	}
	public Gattaca _RunPossibleMutation(double time) {
		StringBuilder MutsObtained = new StringBuilder();
		triNucMuts=new int[triNucProbs.length]; // Must reset to zeroes
		int mutsHappened=0;
		for (int j = 0; j < expectedMuts.length; j++) {
			mutations = PoissonDists[j].nextInt(); // How many mutations for this gene.
			if (mutations>0) {
				mutsHappened+=mutations;
				triNucMuts = new int[triNucProbs.length];
				GattacaMultinomial(triNucProbs, mutations, triNucMuts); // Step One acquire mutation trinucleotide based on 32 mutation types

				for (int tri = 0; tri < triNucMuts.length; tri++) {
					for (int hit = 0; hit < triNucMuts[tri]; hit++) { // Step Two acquire mutation type based on probability of context transition/transversion to get the 96 mutation types.
							theMutOptions[0] = mutTypeProb[(tri * 3) + 0];
							theMutOptions[1] = mutTypeProb[(tri * 3) + 1];
							theMutOptions[2] = mutTypeProb[(tri * 3) + 2];
							NormalizeMutTypeOptions();
							mutTypeHit = new int[3];
							GattacaMultinomial(mutOptionsProbs, 1, mutTypeHit);
							mutIdx = NonZeroIdx(mutTypeHit);
							mutChosen = mutTypes[(tri*3)+mutIdx]; // Which mutation occurs in bases.
							mutPos = basePositions[j][tri][RN.rn.Int(basePositions[j][tri].length)];
							mutChrom = chrom[j];
							mutGene = geneNames[j];
							String mutOut = time + "." + mutGene + "." + mutChosen + "." + mutChrom + "." + mutPos + ";";
							MutsObtained.append(mutOut);
						}
				}
			}

		}

		if(mutsHappened>0){
			String NewPrivateGenome=MutsObtained.toString();
			return(new Gattaca(this, NewPrivateGenome, RN.Double(),RN.Double(),1,RN));
		} else {
			return(this);
		}
	}
private int NonZeroIdx(int[] arr){
		int idx=199;
		int sum=0;
		for (int i = 0; i < arr.length; i++) {
			if(arr[i]!=0){
				idx=i;
			}
			sum+=arr[i];
		}
		if(sum>1 | idx==199){
			System.out.println(arr[0] + "	" + arr[1] + "	" + arr[2]);
			System.out.println(idx);
			throw new IllegalStateException("Multiple mutation types chosen: "+sum);
		}
		return(idx);
	}

	private void BuildPoissons(){
		PoissonDists = new Poisson[expectedMuts.length];
		for (int i = 0; i < expectedMuts.length; i++) {
			Poisson poisson_dist = new Poisson(expectedMuts[i], RNEngine);
			PoissonDists[i] = poisson_dist;
		}
	}

	private static void NormalizeMutTypeOptions(){
		double sum = sumArray(theMutOptions);
		for (int i = 0; i < theMutOptions.length; i++) {
			mutOptionsProbs[i]=theMutOptions[i]/sum;
		}
		double check=1;
		for (int i = 0; i < mutOptionsProbs.length; i++) {
			check-=mutOptionsProbs[i];
		}

		mutOptionsProbs=FixSumToOne(check, mutOptionsProbs);
	}

	private static double sumArray(double[] vals){
		double val=0;
		for (int i = 0; i < vals.length; i++) {
			val+=vals[i];
		}
		return(val);
	}

	private static double[] FixSumToOne(double check, double[] arr){
		if(check!=0){
			for (int i = 0; i < arr.length; i++) {
				if(check<0){
					arr[i]-=check/arr.length;
				} else if(check>0) {
					arr[i]+=check/arr.length;
				}
			}
		}
		return(arr);
	}

	public void GattacaMultinomial(double[] probabilities, int n, int[] ret) {
		double probRemaining = 1;
		for (int i = 0; i < probabilities.length; i++) {
			double prob=probabilities[i];
			if(n==0){
				ret[i]=0;
				return;
			}
			if(probRemaining-prob==0){
				ret[i]=n;
				return;
			}
			int popSelected=picker.Binomial(n,prob/ probRemaining);
			probRemaining -=prob;
			n-=popSelected;
			ret[i]=popSelected;
		}
	}

	// Parses Base Mutation Information
	private static long[][][] ParseTriNucFile(){
		FileIO reader = new FileIO(BaseIndexFile, "r");
		ArrayList<long[]> data = new ArrayList<> (reader.ReadLongs(","));
		long[][][] BaseIndexes = new long[data.size()/32][32][];
		for (int i = 0; i < data.size(); i+=32) {
			BaseIndexes[i/32][0] = data.get(i+0);
			BaseIndexes[i/32][1] = data.get(i+1);
			BaseIndexes[i/32][2] = data.get(i+2);
			BaseIndexes[i/32][3] = data.get(i+3);
			BaseIndexes[i/32][4] = data.get(i+4);
			BaseIndexes[i/32][5] = data.get(i+5);
			BaseIndexes[i/32][6] = data.get(i+6);
			BaseIndexes[i/32][7] = data.get(i+7);
			BaseIndexes[i/32][8] = data.get(i+8);
			BaseIndexes[i/32][9] = data.get(i+9);
			BaseIndexes[i/32][10] = data.get(i+10);
			BaseIndexes[i/32][11] = data.get(i+11);
			BaseIndexes[i/32][12] = data.get(i+12);
			BaseIndexes[i/32][13] = data.get(i+13);
			BaseIndexes[i/32][14] = data.get(i+14);
			BaseIndexes[i/32][15] = data.get(i+15);
			BaseIndexes[i/32][16] = data.get(i+16);
			BaseIndexes[i/32][17] = data.get(i+17);
			BaseIndexes[i/32][18] = data.get(i+18);
			BaseIndexes[i/32][19] = data.get(i+19);
			BaseIndexes[i/32][20] = data.get(i+20);
			BaseIndexes[i/32][21] = data.get(i+21);
			BaseIndexes[i/32][22] = data.get(i+22);
			BaseIndexes[i/32][23] = data.get(i+23);
			BaseIndexes[i/32][24] = data.get(i+24);
			BaseIndexes[i/32][25] = data.get(i+25);
			BaseIndexes[i/32][26] = data.get(i+26);
			BaseIndexes[i/32][27] = data.get(i+27);
			BaseIndexes[i/32][28] = data.get(i+28);
			BaseIndexes[i/32][29] = data.get(i+29);
			BaseIndexes[i/32][30] = data.get(i+30);
			BaseIndexes[i/32][31] = data.get(i+31);
		}
		return BaseIndexes;
	}
}