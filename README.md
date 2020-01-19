# CRISPRScreenProcessing
Function to identify novel sgRNAs in CRISPR organoid screen.
This function has been developed for the analysis the screen of the following study: 

Title: Genome-scale CRISPR screening in organoids identifies synergistic tumor-suppressor activities that lead to TGFß resistance.

Authors: Till Ringel <sup>1</sup> , Nina Frey <sup>1</sup>, Femke Ringnalda <sup>1</sup>, Sharan Janjuha <sup>1</sup>, Sarah Cherkaoui <sup>2</sup>, Stefan Butz <sup>3</sup>, Sumana Srivatsa <sup>4</sup>, Martin Pirkl <sup>4</sup>, Giancarlo Russo <sup>5</sup>, Lukas Villiger <sup>1</sup>, Gerhard Rogler <sup>6</sup>, Hans Clevers <sup>7.8</sup>, Niko Beerenwinkel <sup>4</sup>, Nicola Zamboni <sup>2</sup>, Tuncay Baubec <sup>3</sup> and Gerald Schwank <sup>1</sup>

1. Institute of Molecular Health Sciences, ETH Zurich, Switzerland.
2. Institute of Molecular Systems Biology, ETH Zurich, Switzerland.
3. Department of Molecular Mechanisms of Disease, University of Zurich, Switzerland.
4. Department of Biosystems Science and Engineering, ETH Zurich, Switzerland.
5. Functional Genomics Center Zurich, University of Zurich, ETH Zurich, Switzerland.
6. Department of Gastroenterology and Hepatology, University Hospital Zurich, Switzerland.
7. Oncode Institute, Hubrecht Institute, Royal Netherlands Academy of Arts and Sciences and University Medical Center Utrecht, Utrecht, Netherlands.
8. Princess Máxima Center for Pediatric Oncology, Utrecht, Netherlands.

# Description
Because of the nature of our CRISPR screening approach, in which we do not compare treated and non-treated pools, we decided to use the information of abundant positive controls to identify novel functional sgRNAs. We first selected positive control clones with high read counts for sgRNAs targeting known positive regulators. For these control clones, we identified integrated sgRNAs, removed background sgRNAs reads, and applied the learned pattern to analyze integrated sgRNAs in novel clones. 

Processing function steps:

0. Pre-step. Manual identification of positive control organoid clones with high read counts for sgRNAs targeting known positive regulators. In the TGFß screens, the positive regulators were APC, AXIN, TGFBR1/2
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

Pre-step. Manual identification of positive control organoid clones with high read counts for sgRNAs targeting known positive regulators. In the TGFß screens, the positive regulators were APC, AXIN, TGFBR1/2. An example of such file provided [here](https://github.com/cherkaos/CRISPRScreenProcessing/blob/master/tests/testthat/APK-1-and-2-final.txt).

```
# Running CRISPRScreenProcessing

   inputfile="APK-1-and-2-final.txt"
   controlStart="sample35" 
   controlEnd="sample19-II"
   resultfile <- "APK-1-and-2-final-result.txt"
   maxsgRNA=15
   minReadCount=100
   zscore=TRUE
   screenProcessing(inputfile,controlStart,controlEnd, resultfile, 15,100,TRUE)

```
