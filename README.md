# RegionalAnalysis
![bash](https://img.shields.io/badge/language-bash-green)
![Python](https://img.shields.io/badge/language-Python-blue)
![R](https://img.shields.io/badge/language-R-red)

This repository details scripts and analysis of Chapter Three of my PhD Thesis focusing on the regional diversity of *S. aureus* 
All bioinformatic analysis was conducted on the New Zealand eScience Infrastructure [NeSI](https://github.com/nesi). FastQC and Kraken2 was used for QC of the samples. Nullarbor was used for genome analysis, SpaTyper and AGRvate was used to determine *spa* and *agr* types. 

## FastQC 
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to check for QC of the samples for adaptor content and sequence quality
```#!/bin/bash -e
#SBATCH --cpus-per-task=8 --mem 50Gb --time 1:00:00 -J FASTQC_EV

module load FastQC/0.11.9

FASTQC = /pathtorawreads
OUTPUT_DIR = /processedreadsdirectory/FASTQC/

mkdir -p $OUTPUT_DIR

for file in $FASTQC/*.fq.gz;
do
if [-f $FILE ];
then echo "Running FastQC on $FILE"
fastqc $FILE -o $OUTPUTDIR
else
echo "No FASTQ files found in $FASTQ"
fi
done
echo "FastQC analysis completed for all samples"
```

## Kraken2 
The default Kraken2 installed with [Nullarbor](https://github.com/tseemann/nullarbor) was used for species identification. If an isolate had a *k-mer* match ≤70% to *S. aureus* and/or contained more than 5% *k-mer* matches to another bacterial species these samples will be excluded. 

```#!/bin/bash
CONF=0.5
DATA=`find ${1} -name "*R1*.fq.gz"`
mkdir -p Kraken_Reports/
OUTPUT= Kraken_Reports/

for i in ${DAT[@]}
do
        FQ1=${i}
        FQ2=`echo ${i} | sed s/R1/R2/`
        echo $FQ1
        echo $FQ2
        base=$(basename ${i}  _R1.fq.gz)
        echo "basename is ${base}"
        echo "#!/bin/bash" > tmp.sh
        echo "set -x; module load nullarbor/2.0.20191013; kraken2 --confidence ${CONF} --report ${OUTPUT}/${base}.report --threads 24 --paired ${FQ1} ${FQ2} > ${OUTPUT}/${base}.log" >> tmp.sh
        sbatch -J KrakenIllumina_EV --mem 50G --time 0-1 -c 24 tmp.sh
        sleep 0.5
done
```
## Nullarbor
[Nullarbor](https://github.com/tseemann/nullarbor) is the genomic analysis pipeline chosen for bovine *S. aureus* analysis. In the nullarbor.pl file --mincov was specified as 55 and --minid as 90. The script was run on NESI. 

```#!/bin/bash -e
#SBATCH --cpus-per-task=20 --mem 160Gb --time 166:00:00 -J nullarbor_EV
#Alter text file for purpose - bovine or ST1 analysis
Ref_Bovine = SaureusRF122.fa
##Ref_ST1 = 23EV612.fa

mkdir -p $TMPDIR
module purge
module load nullarbor/2.0.20191013

nullarbor.pl --name StaphAureusIsolates --ref $Ref_Bovine \
--input Bovine_Isolates.txt \
--outdir Bovine_Analysis --force \
--cpus $SLURM_CPUS_PER_TASK --run --mlst saureus --trim --taxoner-opt confidence 0.5 --treebuilder iqtree --treebuilder-opt "-wbtl -st DNA -b 1000 --alrt 1000" --annotator-opt "--fullanno --compliant"
```

## Defining *Spa* Types 
*Spa* Types were determined using [SpaTyper](https://github.com/HCGB-IGTP/spaTyper). The list of genomes is made from the contigs.fa created by Nullarbor. 
```ls *.fa > list_of_genomes.txt ##Creating a list of genomes from all fastq files 
sed -i '1 i\spaTyper  -f ' list_of_genomes.txt
echo "--output spa.xls" &>> list_of_genomes.txt
tr '\n' ' ' < list_of_genomes.txt > spa_command.txt
chmod +x spa_command.txt
./ spa_command.txt
```
## Defining *agr* types
[AGRvate](https://github.com/VishnuRaghuram94/AgrVATE) was used to identify each genomes *agr* type. A conda environment was created then the contig fasta files was used as input for AGRvate. 

## Calculating Statistical Significance Antimicrobial Resistance, Teat Spray Type and Application 
Python script was used to calculate if genes were statistically significant [StatisticallySignficanceGenes.py](https://github.com/emv6/RegionalAnalysis/blob/main/StatisticallySignificant.py) using Chi-Square/Fisher's Exact based on presence/absence data of each detected gene in 309 bovine *S. aureus* isolates 

## Heat Maps for AMR genes per region and per sequence type (ST) 
[R](https://github.com/emv6/RegionalAnalysis/blob/main/Gene_Heatmap.Rmd) script used to generate heatmaps for presence/absence of *blaZ* *fosB* *fusC*, *qacA* and *qacB* across the NZ regions and sequence types (STs). 

## Box Plot generation from SNP matrix of two Taranaki study farms. 
R() script used the pairwise SNP matrix from Nullarbor as input for the the two Taranaki farms (TNK11 and TNK24) to study the genetic diversity between the two farm groups seen on the phylogenetic tree on the phylogenetic tree.  





