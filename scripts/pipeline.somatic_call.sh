#!/bin/bash
source activate somatic

inputdir="$1"
outputdir="$2"

threadN="$3"
jobN="$4"
queueN="$5"
sep="$6"
genome="/mnt/mone/PMI/CH/Reference/Genome/genome.fa"
#targetBed="/mnt/mone/PMI/CH/Reference/On_Target_Bed/Targets-XGEN.66F814ECD1264A978BE2D43A40BDAB94.bed"
preprocess="all"

##Preprocessing
echo "*********** PREPROCESSING(FASTQC, TRIM-GALORE, BWA-MEM, SAMTOOLS) ***********"
soutputdir="$outputdir/01.Alignment/"
samples=['105', '387', '338', '336', '275', '284', '295', '309', '373', '303']

for s in $samples; do
	python2.7 /mnt/mone/PriFiles/sunyme95/scripts/run_preprocess.py -I $inputdir -S $s -O $soutputdir -p $preprocess -t $threadN -j $jobN -q $queueN -s $sep
done

##Run BQSR
echo "*********** Run GATK Base Quality Score Recalibration(BQSR) ***********"
for s in $samples; do
	sinputdir="$soutputdir/dpmrk/$s/"
	soutputdir="$outputdir/01.Alignment/BQSR/$s/"
	echo "python2.7 /mnt/mone/PriFiles/sunyme95/scripts/run_bqsr.py $sinputdir $soutputdir $genome $threadN $jobN $queueN"
	python2.7 /mnt/mone/PriFiles/sunyme95/scripts/run_bqsr.py $sinputdir $soutputdir $genome $threadN $jobN $queueN
done

##Run Mutect2
echo "*********** Run Mutect2 and Filtration steps ***********"
for s in $samples; do
	sinputdir="$soutputdir/BAM/$s/"
	soutputdir="$outputdir/02.Variant_Calling/$s/"
	python2.7 /mnt/mone/PriFiles/sunyme95/scripts/run_mutect2.py -I $sinputdir -O $soutputdir
done

##Run Annotation (Annovar & snpEff)
#echo "*********** Annotate Somatic variants(ANNOVAR) ***********"
#for s in $samples; do
#	sinputdir="$soutputdir/3.filt/$s/"
#	soutputdir="$soutputdir/5.annotation/$s/"
#done

#anno="annovar"
#python2.7 /mnt/mone/PriFiles/sunyme95/scripts/run_annotation.py $sinputdir $soutputdir $targetBed all

 #source activate py2
 #
 ###Filter Somatics
 #sinputdir="$soutputdir/annovar-snpEff/"
 #soutputdir2="$soutputdir/somatics/"
 #python2.7 /mnt/mone/PriFiles/sunyme95/scripts/parse/filtvcf_somatic.py $sinputdir $soutputdir2 $targetBed
 #
 ##Merge VCF files
 #sinputdir="$soutputdir/somatics/"
 #outfile="$soutputdir/somatics.merged.tmp"
 #python2.7 /mnt/mone/PriFiles/sunyme95/scripts/parse/joint_call.py $sinputdir $outfile
 #
 ##sort by chr & position
 #cat $outfile | (sed -u 1q; sort -k1,1V -k2,2n) > $soutputdir/somatics.merged.vcf
 #rm $outfile
 #
 ##Make Matrix
 #infile="$soutputdir/somatics.merged.vcf"
 #matrixfile="$soutputdir/somatics.merged.mtx.txt"
 #python /mnt/mone/PriFiles/sunyme95/scripts/parse/annoVcf2mtx.py $infile $targetBed $matrixfile
 #
 ##Make MAF file for Oncoplot
 #sinputdir="$soutputdir/somatics/"
 #mafdir="$outputdir/MAF/"
 #python /mnt/mone/PriFiles/sunyme95/scripts/run_vcf2maf.py $sinputdir $mafdir
 #
 echo "SOMATIC VARIANT CALLING PIPELINE COMPLETED!!"
