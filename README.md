# CRISPRScreenProcessing
Function to identify novel sgRNAs in CRISPR organoid screen.
This function has been developed for the analysis the screen of the following study: 

Title: Genome-scale CRISPR screening in organoids identifies synergistic tumor-suppressor activities that lead to TGFß resistance.

Authors: Till Ringel <sup>1,2</sup>, Nina Frey <sup>1,2</sup>, Femke Ringnalda <sup>1,2</sup>, Sharan Janjuha <sup>1,2</sup>, Sarah Cherkaoui <sup>4</sup>, Stefan Butz <sup>5</sup>, Sumana Srivatsa <sup>6</sup>, Martin Pirkl <sup>6</sup>, Giancarlo Russo <sup>7</sup>, Lukas Villiger <sup>1,2</sup>, Gerhard Rogler <sup>8</sup>, Hans Clevers <sup>3,9</sup>, Niko Beerenwinkel <sup>6</sup>, Nicola Zamboni <sup>4</sup>, Tuncay Baubec <sup>5</sup> and Gerald Schwank <sup>1,2,10,*</sup>

<sup>1</sup>Institute of Molecular Health Sciences, ETH Zurich, Switzerland.

<sup>2</sup>Department of Pharmacology and Toxicology, University of Zurich, Switzerland.

<sup>3</sup>University Medical Center (UMC) Utrecht, Utrecht, Netherlands.

<sup>4</sup>Institute of Molecular Systems Biology, ETH Zurich, Switzerland.

<sup>5</sup>Department of Molecular Mechanisms of Disease, University of Zurich, Switzerland.

<sup>6</sup>Department of Biosystems Science and Engineering, ETH Zurich, Switzerland.

<sup>7</sup>Functional Genomics Center Zurich, University of Zurich, ETH Zurich, Switzerland.

<sup>8</sup>Department of Gastroenterology and Hepatology, University Hospital Zurich, Switzerland.

<sup>9</sup>Hubrecht Institute, Royal Netherlands Academy of Arts and Sciences (KNAW), Utrecht, Netherlands

<sup>10</sup>Lead Contact

*Correspondence: [schwank@pharma.uzh.ch](mailto:schwank@pharma.uzh.ch)

# Description
Because of the nature of our CRISPR screening approach, in which we do not compare treated and non-treated pools, we decided to use the information of abundant positive controls to identify novel functional sgRNAs. We first selected positive control clones with high read counts for sgRNAs targeting known positive regulators. For these control clones, we identified integrated sgRNAs, removed background sgRNAs reads, and applied the learned pattern to analyze integrated sgRNAs in novel clones. 

Processing function steps:

0. Pre-step. Manual identification of positive control organoid clones with high read counts for sgRNAs targeting known positive regulators.
1. For all integrations, Z-scores of the reads are calculated to make read count values comparable to each other.
2. For each control clone, the processing function ranks sgRNAs from high to low and identifies largest drop in read count (fold change) between two consecutive sgRNAs.
3. The smallest read count drop (fold change) among control clones is selected as minimum threshold for integrated sgRNAs.
4. Identification of new sgRNAs in non-control clones: only clones with sgRNAs read count drop bigger than the threshold are selected by processing function.
5. Generation of summary table containing control clones and non-control clones (new candidates) with identified integrated sgRNAs. 


# Installation
```
# If you do not have devtools installed
install.packages("devtools")

# Install CRISPRScreenProcessing 
library(devtools)
devtools::install_github("cherkaos/CRISPRScreenProcessing")
library("CRISPRScreenProcessing")
```
# Running the function for the identification of sgRNAs in TGFß screens

Pre-step. Manual identification of positive control organoid clones with high read counts for sgRNAs targeting known positive regulators. In the TGFß screens, the positive regulators were APC, AXIN, TGFBR1/2. An example of such file provided [here](https://github.com/cherkaos/CRISPRScreenProcessing/blob/master/tests/testthat/APK-1-and-2-final.txt). The file should be a txt or csv file (comma-separated) in your working directory.

```
# Running CRISPRScreenProcessing
   
   inputfile="APK-1-and-2-final.txt" 
   controlStart="sample35" 
   controlEnd="sample19-II"
   resultfile <- "APK-1-and-2-final-result.txt"
   maxsgRNA=15
   minReadCount=100
   zscore=TRUE
   orderOutput=TRUE
   shortOutput=TRUE
   screenProcessing(inputfile,controlStart,controlEnd,resultfile,maxsgRNA,minReadCount,zscore,orderOutput,shortOutput)

```
# Output

The output is a file containing control clones and succesfully select novel clones. Each clones has 3 columns, where you can find on the information of the integrated sgRNAs and its read count. An example of such output can be found [here](https://github.com/cherkaos/CRISPRScreenProcessing/blob/master/tests/testthat/APK-1-and-2-final.txt).


