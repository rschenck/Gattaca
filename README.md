# Gattaca

Gattaca is a method for tracking base pair resolution data within agent based simulations. It consists of three parts the 
1. [Setup](#part-1-setup)
2. [Simulations](#part-2-simulations)
3. [Analysis](#part-3-analysis)

## Pre-requisites

Gattaca depends on having a proper reference genome. That genome must match a snpEff database. The following are dependencies.

1. [snpEff](http://snpeff.sourceforge.net)
2. Gunzip compressed primary genome assembly (recommended GRCh37 or GRCh38 from [Ensembl](https://www.ensembl.org/index.html), located from ftp download)
3. [biopython](https://biopython.org/wiki/Download)
4. [HAL](https://halloworld.org)
5. [Colt jar](https://dst.lbl.gov/ACSSoftware/colt/) for use within HAL.

## Part 1: Setup

Change the directory locations for the pre-requisites in the ```usr_path.ini``` file. This will be the location of snpEff and the reference genome.

An example execution from within the Gattaca directory:

```angular2html
python Gattaca.py --geneList ./tests/TestGenes.txt --genome=GRCh37.75 --contextFile=./Tests/MutContext.txt --mutRate=3.2E-9
```

#### Part 1 Output

Gattaca will yield two files that must be placed in the executable directory (generally /src) for simulations within HAL.
1. Gattaca.java
2. triNucsPos.csv

Once these are placed within your simulation project move to Part 2.

## Part 2: Simulations

While other java frameworks can be used or a port of the java code could be constructed for other simulation frameworks Gattaca is designed to work flawlessly within HAL. Thus an additional pre-requisite is [HAL](https://halloworld.org). Gattaca uses some random number generators and distributions from the Colt library ([Colt jar](https://dst.lbl.gov/ACSSoftware/colt/)) as well, which can be placed in the HAL lib directory. This can easily be linked from within ideaJ (instructions on linking libraries can be found from [HAL](https://halloworld.org) or ideaJ)

Integration into HAL requires a few steps:
1. Place the two output files from part 1 within the scope of your main executible class.
2. From within your main function initialize the first clone, this provides a root clone for downstream tracking (also speeds up your simulations):
```angular2html
// The Gattaca constructor requires: parent, String, Hue, Saturation, Value, Rand RNG
Gattaca initClone = new Gattaca(null, "", 1, 0, 0.3, RN);
```
3: Pass initClone into your Grid class. You will initialize clone1 using this, so that your first cell is constructed from within your grid using:
```angular2html
Gattaca clone1=new Gattaca(this.clone0, "",1,0,0.3, RN);
Cell c=NewAgentSQ(xpt, ypt).Init(clone1, BIRTHPROBABILITY); // Birth cell with genome clone1
c.genome.IncPop(); // Increases genome population by 1
```
4: Make sure that your cell class has a Gattaca variable.
5: Increase and decrease the population of the genome using a cells genome DecPop() and IncPop() commands.
6: From within your main model step. You can choose to record clones at any timepoint:
```angular2html
initClone.RecordClones(G.GetTick());
```
7: Output information to file:
```angular2html
String[] AttributesList = new String[]{"Genome", "H", "S", "V"};
initClone.OutputClonesToCSV("/Users/rschenck/Dropbox/GATTACA/Gattaca/tests/GattacaEx/gattaca_output_fullyseeded." + Integer.toString(CON.SEED) + ".csv", AttributesList, (Gattaca g) -> {
    return GetAttributes(g);
}, 0);

// Function to retrieve the attributes of your choice.
public static String GetAttributes(Gattaca root) {
        return root.PrivateGenome + "," + Double.toString(root.h) + "," + Double.toString(root.s)+ "," + Double.toString(root.v);
}
```

## Part 3: Analysis

Ideally, simulations will be ran in replicate for downstream statistical analysis, but a single simulation can also be handled by Gattaca analysis.

```angular2html
python value
```

## Plan

1d, 2d, 3d, gland vs stratified epithelia.