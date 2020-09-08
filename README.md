# Gattaca

Gattaca is a method for tracking base pair resolution data within agent based simulations. It consists of three parts the setup, the execution, and the analysis.

## Pre-requisites

Gattaca depends on having a proper reference genome. That genome must match a snpEff database. The following are dependencies.

1. [snpEff](http://snpeff.sourceforge.net)
2. Gunzip compressed primary genome assembly (recommended GRCh37 or GRCh38 from [Ensembl](https://www.ensembl.org/index.html))

## Part 1: Setup

An example execution from within the Gattaca directory:

```angular2html
python Gattaca.py --geneList ./tests/TestGenes.txt --genome=GRCh37.75 --contextFile=./Tests/MutContext.txt --mutRate=3.2E-9
```

#### Part 1 Output

Gattaca will yield two files that must be placed in the executable directory (generally /src) for simulations within HAL.
1. Gattaca.java
2. triNucsPos.csv

Once these are placed within your simulation project move to Part 2.

## Part 2: Execution

While other java frameworks can be used or a port of the java code could be constructed for other simulation frameworks Gattaca is designed to work flawlessly within HAL. Thus an additional pre-requisite is [HAL](https://halloworld.org). Gattaca uses some random number generators and distributions from the Colt library ([Colt jar](https://dst.lbl.gov/ACSSoftware/colt/)) as well, which can be placed in the HAL lib directory. This can easily be linked from within ideaJ (instructions on linking libraries can be found from [HAL](https://halloworld.org) or ideaJ)



## Plan

1d, 2d, 3d, gland vs stratified epithelia.