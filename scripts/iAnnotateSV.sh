#!/bin/bash

#MUST BE TUMOR NORMAL PAIRED CALLING
#BASEDIR=$(dirname "$0")
BASEDIR=$(readlink -f "$0")
BASEDIR=`echo $BASEDIR | sed 's,/iAnnotateSV.sh,,g'`

#'somatic.vcf.gz from manta output'
vcf_gz=$1
# 'sample name'
prefix=$2
# 'output directory'
outdir=$3
# manta directory 
mantadir=$4
# fasta reference
fasta=$5

echo $BASEDIR
gitdir=`echo $BASEDIR | sed 's,/scripts$,,g'`
# 'resource files for iannotatesv'
resourcedir="${gitdir}/data"
# iannotateSV directory /ifs/work/bergerm1/zhengy1/RET_all/IannotateSV_testing/iAnnotateSV
iasvdir="${gitdir}/iAnnotateSV"

# get rid of end '/' in directories
outdir=`echo $outdir | sed -r 's,/$,,g'`
resourcedir=`echo $resourcedir | sed -r 's,/$,,g'`

vcf_manta="$outdir/`echo $vcf_gz | sed -r 's,.*./,,g' | sed -r 's,.vcf.gz,.vcf,g'`"
echo 'unpacking vcf gz file from manta'
zcat $vcf_gz > $vcf_manta

vcf_input=${vcf_manta/.vcf/_inv_corrected.vcf}
echo 'correcting inversions in manta vcf'
${mantadir}/libexec/convertInversion.py `which samtools` \
  $fasta   $vcf_manta > $vcf_input

edited_vcf="$outdir/`echo $vcf_input | sed -r 's,.*./,,g' | sed -r 's,.vcf,_edited.vcf,g'`"
echo 'Adding connection type and making changes to manta VCF'
Rscript "${BASEDIR}/vcf_ct_edits.R" \
  -v $vcf_input -o $edited_vcf

echo 'Converting manta VCF into tab file for iAnnotateSV'
python "${BASEDIR}/manta_vcf2tab.py" \
	-i $edited_vcf -o $outdir

tab_txt="$outdir/`echo $edited_vcf | sed -r 's,.*./,,g' | sed -r 's,.vcf,.tab,g'`"
echo $tab_txt
echo 'Running iAnnotateSV'
# python "${BASEDIR}/iAnnotateSV.py" \
#python /home/shahr2/git/iAnnotateSV/iAnnotateSV/iAnnotateSV.py \
python "${iasvdir}/iAnnotateSV/iAnnotateSV.py" \
	-r hg19 -d 3000 -i $tab_txt -ofp $prefix -o $outdir \
	-c "${resourcedir}/canonical_transcripts_cv5.txt" \
	-u "${resourcedir}/hg19.uniprot.spAnnot.table.txt" \
	-cc "${resourcedir}/cancer_gene_census.tsv" \
	-cct "${resourcedir}/cosmic_fusion_counts.tsv" \
	-rr "${resourcedir}/hg19_repeatRegion.tsv" \
	-dgv "${resourcedir}/hg19_DGv_Annotation.tsv"

annot_txt="$outdir/`echo $prefix`_Annotated.txt"
annot_output=`echo $annot_txt | sed -r 's,_Annotated.txt,_Annotated_Evidence.txt,g'`

echo 'Merging vcf and iAnnotateSV output'
Rscript "${BASEDIR}/merge_vcf_tab.R" -t $annot_txt \
	-v $edited_vcf -o $annot_output

echo 'Adding DMP IMPACT fusion frequency'
java -server -Xms4g -Xmx4g -cp /home/patelju1/software/BioinfoUtils-1.0.0.jar \
	org.mskcc.juber.commands.AnnotateSVResults \
	$annot_output \
	"${resourcedir}/dmp-intragenic-white-list.txt" \
	"${resourcedir}/dmp-intergenic-white-list.txt"


# from config file	
#ANNOSV:/home/shahr2/git/iAnnotateSV/iAnnotateSV/iAnnotateSV.py
#GENOMEBUILD:hg19
#DISTANCE:3000
#CANONICALTRANSCRIPTFILE:/home/patelju1/resources/icallsv/data/canonical_transcripts_cv5.txt
#UNIPROTFILE:/home/patelju1/resources/icallsv/data/hg19.uniprot.spAnnot.table.txt
#CosmicCensus:/home/patelju1/resources/icallsv/data/cancer_gene_census.tsv
#CosmicFusionCounts:/home/patelju1/resources/icallsv/data/cosmic_fusion_counts.tsv
#RepeatRegionAnnotation:/home/patelju1/resources/icallsv/data/hg19_repeatRegion.tsv
#DGvAnnotations:/home/patelju1/resources/icallsv/data/hg19_DGv_Annotation.tsv
