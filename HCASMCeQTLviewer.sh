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
wget https://www.dropbox.com/s/edm0ykexjmue5yf/reverse.zip?dl=0
gunzip reverse.zip
fi

cd ~/HCASMC

#check if folder with genotypes exist, if not, create

if [ ! -d $GENOTYPES ]
then
mkdir ~/HCASMC/HCASMC_genotypes
fi

#check if subfolder with genotypes exists, if not, download from Dropbox link

if [ ! -f $REV ]
then
wget 
gunzip 
fi
