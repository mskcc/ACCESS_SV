#!/bin/bash

# ----------------------------------------------------------------------------------------------------
# Title:       iAnnotateSV.sh
# Description: Bash wrapper to post-process Manta output and annotate the results with iAnnotateSV
#              - Currently only Tumor - Normal paired calling is supported
# ----------------------------------------------------------------------------------------------------

# enable bash strict mode
set -euo pipefail

# Set default path
BASEDIR=$(readlink -f "$0" | sed 's,/iAnnotateSV.sh,,g')
SAMTOOLS=${$(command -v samtools):-}
PYTHON=${$(command -v python):-}
RSCRIPT=${$(command -v Rscript):-}

# checks
if [ -z "$SAMTOOLS" ] || [ -z "$PYTHON" ]; then
    echo "One or more of samtools or python not defined in path."
    exit 1
fi

# Assign required arguments to variables
# overriding PYTHON with user provided path
read -r vcf_gz prefix outdir mantadir fasta PYTHON <<< "{@:1:6}"
if [ -z "$vcf_gz" ] || [ -z "$prefix" ] || [ -z "$outdir" ] || \
	[ -z "$mantadir" ] || [ -z "$fasta" ] || [ -z "$PYTHON" ]; then
    echo "One or more of the required arguments is missing."
    exit 1
fi

# Declare required paths and remove trailing forward slash
gitdir=$(echo $BASEDIR | sed 's,/scripts$,,g')
resourcedir="${gitdir}/data"
iasvdir="${gitdir}/iAnnotateSV"
outdir=$(echo $outdir | sed -r 's,/$,,g')
resourcedir=$(echo $resourcedir | sed -r 's,/$,,g')

vcf_manta="$outdir/${prefix}_$(echo $vcf_gz | sed -r 's,.*./,,g' | sed -r 's,.vcf.gz,.vcf,g')"
echo 'unpacking vcf gz file from manta'
zcat $vcf_gz > $vcf_manta

vcf_input=${vcf_manta/.vcf/_inv_corrected.vcf}
echo 'correcting inversions in manta vcf'
${mantadir}/libexec/convertInversion.py $SAMTOOLS \
  $fasta   $vcf_manta > $vcf_input

edited_vcf="$outdir/$(echo $vcf_input | sed -r 's,.*./,,g' | sed -r 's,.vcf,_edited.vcf,g')"
echo 'Adding connection type and making changes to manta VCF'
$RSCRIPT "${BASEDIR}/vcf_ct_edits.R" \
  -v $vcf_input -o $edited_vcf

echo 'Converting manta VCF into tab file for iAnnotateSV'
$PYTHON "${BASEDIR}/manta_vcf2tab.py" \
	-i $edited_vcf -o $outdir

tab_txt="$outdir/$(echo $edited_vcf | sed -r 's,.*./,,g' | sed -r 's,.vcf,.tab,g')"
echo 'Running iAnnotateSV'
$PYTHON "${iasvdir}/iAnnotateSV/iAnnotateSV.py" \
	-r hg19 -d 3000 -i $tab_txt -ofp $prefix -o $outdir \
	-c "${resourcedir}/canonical_transcripts_cv5.txt" \
	-u "${resourcedir}/hg19.uniprot.spAnnot.table.txt" \
	-cc "${resourcedir}/cancer_gene_census.tsv" \
	-cct "${resourcedir}/cosmic_fusion_counts.tsv" \
	-rr "${resourcedir}/hg19_repeatRegion.tsv" \
	-dgv "${resourcedir}/hg19_DGv_Annotation.tsv"

annot_txt="$outdir/$(echo $prefix)_Annotated.txt"
annot_output=$(echo $annot_txt | sed -r 's,_Annotated.txt,_Annotated_Evidence.txt,g')

echo 'Merging vcf and iAnnotateSV output'
$RSCRIPT "${BASEDIR}/merge_vcf_tab.R" -t $annot_txt \
	-v $edited_vcf -o $annot_output

echo 'Annotating with DMP IMPACT occurrence count and SplitRead AF'
$RSCRIPT "${BASEDIR}/sv_utilities.R" \
	-c "${resourcedir}/dmp_IMPACT_sv_occurrence_count_summary.txt" \
	-a $annot_output
