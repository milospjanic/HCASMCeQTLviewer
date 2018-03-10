# HCASMCeQTLviewer

HCASMCeQTLviewer is a combined bash/R script to view eQTL box/dot plots for a specific gene/SNP eQTL assocoation in human coronary artery smooth muscle cells, one of the crucial cell types that are involved in atherosclerotic process of the blood vessel.

HCASMCeQTLviewer can be used to quicky check and graph the **directionality of SNP-gene association** in human coronary artery smooth muscle cells.

RNAseq expression and whole genome sequencing data come from our collection of 52 HCASMC cell lines and are part of the upcoming publication. Data repositories used in this tool are private and are available upon request and will be made public upon publication.

# Usage

To run the script download the .sh file

<pre>
wget https://raw.githubusercontent.com/milospjanic/HCASMCeQTLviewer/master/HCASMCeQTLviewer.sh
chmod 755 HCASMCeQTLviewer.sh
</pre>

Place the script in your home or any other folder. The script will create ~/HCASMC as its working folder, and two subfolders: **HCASMC_expr** and **HCASMC_genotypes**. HCASMC_expr will contain per-gene RNAseq read counts for each HCASMC sample, while HCASMC_genotypes contains whole genome sequencing vcf files of HCASMC samples. Script will create subfolder ~/HCASMC/GENE-SNP for each GENE-SNP input combination and place the textual and pdf output files in that folder. Text file **FINAL.txt.cut** contains genotypes and expression levels:

<pre>
29.5578 CC
34.5678 CG
42.2207 CC
38.6143 CG
33.2066 CG
...
</pre>

Pdf outputs are combind box/dot plots with rsID and rsID genotypes on x-axis and GENE normalized expression levels on y-axis (reads per million mapped).

To run the script create the file named **genes** containing the GENE name, and provide it as first argument. Second and third arguments should be SNP and chromosome on which the SNP is located.

<pre>
cat gene
SLC22A4

./HCASMCeQTLviewer.sh gene rs1537373 9 
</pre>

Script will check if all the expression and genotype data sets are present in **HCASMC_expr** and **HCASMC_genotypes** subfolders and if not it will download them from the repository.

# Examples
Example of a strong eQTL for FES gene.

![alt text](https://github.com/milospjanic/HCASMCeQTLviewer/blob/master/rs2521501.png)

Example of an eQTL for SLC22A4 gene.

![alt text](https://github.com/milospjanic/HCASMCeQTLviewer/blob/master/rs273909.png)

Example of a weak association for CDKN2B gene.

![alt text](https://github.com/milospjanic/HCASMCeQTLviewer/blob/master/rs1537373.png)

Example of an eQTL for TCF21 genes.

![alt text](https://github.com/milospjanic/HCASMCeQTLviewer/blob/master/rs12190287.png)
