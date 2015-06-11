[![DOI](https://zenodo.org/badge/13132/darencard/RADpipe.svg)](http://dx.doi.org/10.5281/zenodo.17809)

RADpipe
=======

## Overview:
This repository contains a set of scripts that will pipeline data processing of RADseq data to produce files necessary for analysis. Details of this pipeline and of some supplementary scripts are included below. Please note that some up-front data processing steps are specific to the type of library prep and reads used in my work and may not work properly in instances where alternative library preparations or adapter designs are used. These scripts are designed to process double digest RADseq data with combinatorial barcodes, like that presented in Peterson et al. (2012). These data also include an 8bp unique molecular identifier (UMI) sequence at the beginning of both the forward and reverse (if applicable) reads that enable PCR clones to be excluded. With this in mind, considerable flexibility is included in the pipeline scripts, which may allow other library types to be accommodated. Moreover, anyone with working knowledge of Python can adapt these scripts to their design. I apologize in advance for the inefficient nature of many portions of the code and I know that it would be largely frowned upon by a more seasoned coder, but ultimately it does work (at least for me).

If you decide to use or adapt these scripts for your own work, please be sure to acknowledge them, and all dependencies, in any resulting publications. Full citations have been included below to facilitate this.

Please note that I do not have the time or expertise to provide support, so these scripts are provided as-is, with no guarantee of proper functioning or desirable results. Use with caution!

## Dependencies:
1. Python 2 (tested using v.2.7.3)
2. Stacks up to v.1.19 (tested using v.1.19)
3. BWA (tested using v.0.7.9)
4. SAMtools version 0 (tested using v.0.1.19)
5. BCFtools version 0 (tested using v.0.1.19)
6. VCFtools (tested using v.0.1.12b)
7. R with Base, Utils, Stats, MASS, and RColorBrewer packages installed

## Core Pipeline:
1. process_rawreads.py: Filters PCR clones, trims away 8bp UMI, parses reads for each sample, and quality trims.
2. read_mapping.py: Maps parsed reads from each sample to a specified reference genome and processes mapping files.
3. variant_calling_from_BAM.py: Uses mapping BAM files to call variants and produce a VCF output.

## Other Scripts:
1. sigThreshold_bootstrap.py: Returns a threshold for significance based on bootstrap resampling of a given column (i.e., a population genetic statistic).
2. genotype_from_VCF.py: Returns either genotype likelihood matrices (with format designated by user) or variant alignments for downstream programs.
3. entropyStart.R: Produces MCMC starting points for the Entropy program (Gompert et al. 2014) using output from genotype_from_VCF.py.
4. meta_sort_NGSadmix.py: Formats admixture proportion output from NGSadmix () so it can be manipulated and plotted using admixturePlot.R. Will likely adjust so that alternate outputs can be parsed.
5. admixturePlot.R: Produces admixture bar plot (i.e., "Structure" plot) from formated output from meta_sort_NGSadmix.py for visualization.

## Running the Pipeline:
Given that each script contains detailed usage information, no further details will be provided here for now. I hope to start filling in examples as time permits.

## Acknowledgements:
This pipeline benefited from many discussions with other members of the Castoe lab at the University of Texas at Arlington. Chris Nice (Texas State University) provided guidance on SNP calling and formating VCF files for downstream analyses.

## Citations
##### This Repository:
Card. 2015. RADpipe. GitHub repository. [doi: 10.5281/zenodo.17809](http://doi.org/10.5281/zenodo.17809).

##### ddRADseq:
Peterson, et al. 2012. Double Digest RADseq: An inexpensive method for de novo SNP discovery and genotyping in model and non-model species. PLoS ONE 7 (5): e37135. [doi: 10.1371/journal.pone.0037135](http://doi.org/10.1371/journal.pone.0037135).

##### Stacks:
Catchen, et al. 2011. Stacks: building and genotyping loci de novo from short-read sequences. G3: Genes, Genomes, Genetics 1 (3): 171-182. [doi: 10.1534/g3.111.000240](http://doi.org/10.1534/g3.111.000240).
Catchen, et al. 2013. Stacks: an analysis tool set for population genomics. Molecular Ecology 22 (11): 3124-3140. [doi: 10.1111/mec.12354](http://doi.org/10.1111/mec.12354).

##### BWA/SAMtools/BCFtools:
Li & Durbin. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25 (14): 1754-1760. [doi: 10.1093/bioinformatics/btp324](http://doi.org/10.1093/bioinformatics/btp324).
Li, et al. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25 (16): 2078-2079. [doi: 10.1093/bioinformatics/btp352](http://doi.org/10.1093/bioinformatics/btp352).
Li. 2011. A statistical framework for SNP calling, mutation discovery, association mapping, and population genetical parameter estimation from sequence data. Bioinformatics 27 (21): 2987-2993. [doi: 10.1093/bioinformatics/btr509](http://doi.org/10.1093/bioinformatics/btr509).

##### VCFtools:
Danecek, et al. 2011. The variant call format and VCFtools. Bioinformatics 27 (15): 2156-2158. [doi: 10.1093/bioinformatics/btr330](http://doi.org/10.1093/bioinformatics/btr330).

##### Admixture Analyses:
Gompert, et al. 2014. Admixture and the organization of genetic diversity in a butterfly species complex revealed through common and rare genetic variants. Molecular Ecology 23 (18): 4555-4573. [doi: 10.1111/mec.12811](http://doi.org/10.1111/mec.12811).
Skotte, et al. 2013. Estimating individual admixture proportions from next generation sequencing data. Genetics 195 (3): 693-702. [doi: 10.1534/genetics.113.154138](http://doi.org/10.1534/genetics.113.154138).

##### R:
R Development Core Team. 2015. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. [http://www.r-project.org/](http://www.r-project.org/).
