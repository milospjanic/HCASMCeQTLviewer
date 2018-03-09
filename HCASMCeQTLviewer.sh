
#!/bin/bash

GENE=$(pwd)/$1
EXPR=~/HCASMC/HCASMC_expr
REV=~/HCASMC/HCASMC_expr/reverse
GENOTYPES=~/HCASMC/HCASMC_genotypes/
VCF=~/HCASMC/HCASMC_genotypes/vcf
SNP=$1
CHR=$2

cd 
mkdir ~/HCASMC

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
echo "Unpacking expression files..."
unzip reverse.zip
fi

cd ~/HCASMC

#check if folder with genotypes exist, if not, create

if [ ! -d $GENOTYPES ]
then
mkdir ~/HCASMC/HCASMC_genotypes
fi

#check if subfolder with genotypes exists, if not, download from Dropbox link

cd ~/HCASMC/HCASMC_genotypes/

if [ ! -f $VCF ]
then 
mkdir ~/HCASMC/HCASMC_genotypes/vcf
cd ~/HCASMC/HCASMC_genotypes/vcf
wget https://www.dropbox.com/s/nnytxlbx1v0gh8y/phased_and_imputed.tar
tar -zvf phased_and_imputed.tar
echo "Unpacking genome vcf files..."
gunzip phased_and_imputed*
fi

cd ~/HCASMC/HCASMC_expr

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

awk 'NR==FNR {h[$2] = $1; h2[$2] = $2; next} {print h[$1]}' id_merge.txt $GENE > genename

#remove temporary files

rm id_merge.txt
rm script.r

#get gene counts for gene of interest

while read line; do
                set $line
                find . -name *gene.count | xargs grep $line > COUNTS.txt
done < genename

sed -E "s/^.\/reverse\///g" COUNTS.txt | sed -E "s/\/.*\.[0-9]*//g" > COUNTS.txt.cut

#get total gene counts per sample

find . -name *gene.count | xargs -I % awk 'BEGIN {FS = " "} ; {sum+=$2} END {print sum}' % > TOTAL.txt
find . -name *gene.count | xargs -I % echo % > SAMPLES.txt 

sed -E "s/^.\/reverse\///g" SAMPLES.txt | sed -E "s/\/.*\.[0-9a-zA-Z]*//g" > SAMPLES.txt.cut

paste SAMPLES.txt.cut TOTAL.txt > TOTALCOUNTS.txt

awk 'NR==FNR {h[$1] = $0; next} {if(h[$1]) print h[$1]"\t"$0}' COUNTS.txt.cut TOTALCOUNTS.txt > TABLE.txt

awk '{print $1 "\t" $2/$4*1000000}' TABLE.txt > TABLE.RPM.txt

#create genotype files

grep $SNP phased_and_imputed.chr$CHR.vcf > SNP.txt
grep CHROM phased_and_imputed.chr$CHR.vcf > HEADER.txt


#clean genotype files
sed -i -E 's/:[0-9.]*:[0-9.]*,[0-9.]*,[0-9.]*//g' SNP.txt 

#greb reference and alternative aleles
REF="$(awk '{printf $4}' SNP.txt)"
ALT="$(awk '{printf $5}' SNP.txt)"

sed -i "s/0|0/$REF$REF/g" SNP.txt
sed -i "s/0|1/$REF$ALT/g" SNP.txt
sed -i "s/1|1/$ALT$ALT/g" SNP.txt

cat HEADER.txt SNP.txt > GENOTYPES.txt

#write script for trasposing a matrix in awk

echo "{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = \$i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str\" \"a[i,j];
        }
        print str
    }
}
" > transpose.awk 

sed -i 's/#CHROM/CHROM/g' GENOTYPES.txt
awk -f transpose.awk GENOTYPES.txt > GENOTYPES.txt.tr

#compile two tables 

awk 'NR==FNR {h[$1] = $0; next} {if(h[$1]) print h[$1]"\t"$0}' TABLE.RPM.txt GENOTYPES.txt.tr

#Rcode
#library(ggplot2)
#data<-read.table (file="rs12190287", sep="\t",head=T)
#colnames(data)<-c("Genotype","TCF21 expression")
#p <- ggplot(data, aes(x=data$"Genotype",y=data$"TCF21 expression")) + geom_boxplot()
#pdf("rs12190287.pdf")
#p+scale_x_discrete(limits=c("GG", "CG", "CC"))+geom_jitter(shape=16, position=position_jitter(0.2))+theme(axis.text=element_text(size=24),axis.title=element_text(size=26))+ labs(title = "rs2327433", x="Genotype", y="TCF21 expression") + theme(plot.title = element_text(size = rel(2)))
#dev.off()


