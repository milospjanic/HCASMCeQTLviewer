
#!/bin/bash

EXPR=~/HCASMC/HCASMC_expr
REV=~/HCASMC/HCASMC_expr/reverse
GENOTYPES=~/HCASMC/HCASMC_genotypes/


#check if folder with expression values exist, if not, create

if [ ! -d $EXPR ]
then
mkdir ~/HCASMC/HCASMC_expr
fi

cd ~/HCASMC/HCASMC_expr

#check if reverse folder with expression levels counted on reverse strand exists, if not, download from Dropbox link

if [ ! -f $REV ]
then
wget https://www.dropbox.com/s/edm0ykexjmue5yf/reverse.zip
gunzip reverse.zip
fi

cd ~/HCASMC

#check if folder with genotypes exist, if not, create

if [ ! -d $GENOTYPES ]
then
mkdir ~/HCASMC/HCASMC_genotypes
fi

#check if subfolder with genotypes exists, if not, download from Dropbox link

#if [ ! -f $REV ]
#then
#wget 
#gunzip 
#fi


#write R script to get ENSEMBL id, needs biomaRt in R

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")
ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"hsapiens_gene_ensembl\", host=\"grch37.ensembl.org\")
id_merge = getBM(attributes=c(\"ensembl_gene_id\",\"external_gene_name\"),mart=ensembl)
write.table(id_merge, file=\"id_merge.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)
" > script.r

#run R script

chmod 775 script.r
./script.r

#Use awk to append gene names

awk 'NR==FNR {h[$1] = $1; h2[$1] = $2; next} {print h[$1]"\n"}' id_merge.txt $1 > $1.genename

#remove temporary files

rm id_merge.txt
rm script.r

#get gene counts

cd ~/HCASMC/HCASMC_expr
while read line; do
                set $line
                find . -name *gene.count | xargs grep $line > COUNTS.txt
                sed -E "s/^.\///g" COUNTS.txt | sed -E "s/\/.*\.[0-9]*//g" 
done < $1.genename


