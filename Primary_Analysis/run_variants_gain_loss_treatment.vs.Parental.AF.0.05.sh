#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"

data_path=DATA_PATH
project=PROJECT
subject=SUBJECT
sample_ctl="PARENTAL_SAMPLE"
sample_treatment="TREATED_SAMPLE"

vcf_ctl=${data_path}/${project}/${subject}/${sample_ctl}/${sample_ctl}.aligned.duplicates_marked.recalibrated.targeted_sequencing.somvarius.tumor_only.sorted.rm_dbsnps.vcf

# index the vcf if not already
#if [ -e ${vcf_ctl}.idx ];then 
/data/ldiaz/luol2/gatk-4.1.2.0/gatk IndexFeatureFile -F ${vcf_ctl}
#fi

vcf_treatment=${data_path}/${project}/${subject}/${sample_treatment}/${sample_treatment}.aligned.duplicates_marked.recalibrated.targeted_sequencing.somvarius.tumor_only.sorted.rm_dbsnps.vcf

# index the vcf if not already
#if [ -e ${vcf_treatment}.idx ];then 
/data/ldiaz/luol2/gatk-4.1.2.0/gatk IndexFeatureFile -F ${vcf_treatment}
#fi

#######################################################################
# Variant gain against Parental
####################################################################### 

# grab difference as treatment variant gain from Parental
command_mem="1000"
/home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${vcf_treatment} --discordance ${vcf_ctl} -O ${data_path}/${project}/${subject}/${sample_ctl}_vs_${sample_treatment}.vcf
# rename the vcf to Gain.xxxxxxx
tab=`echo -ne "\t"`;sed -i "s/\(CHROM.*FORMAT\).*/\1${tab}Gain.${sample_treatment}/" ${data_path}/${project}/${subject}/${sample_ctl}_vs_${sample_treatment}.vcf

module load singularity/3.1.1
VEPInputVcf=${sample_ctl}_vs_${sample_treatment}.vcf
VEPOutputVcf=${sample_ctl}_vs_${sample_treatment}.VEP.ann.vcf
VEPOutputStat=${sample_ctl}_vs_${sample_treatment}.VEP.ann.summary.html

#Perform VEP using local cache (Homo_Sapiens)
singularity run --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
vep -i /data/${VEPInputVcf} \
    -o /data/${VEPOutputVcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --stats_file /data/${VEPOutputStat}


#Perform VEP summarization after manual filtration
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' ${data_path}/${project}/${subject}/${VEPInputVcf} > ${data_path}/${project}/${subject}/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.vcf} 
singularity exec --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.vcf} \
    -o /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.VEP.ann.vcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --stats_file /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.VEP.ann.summary.html}

# VEP output in tab format
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.vcf} \
    -o /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.VEP.ann.txt} \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Gene,CANONICAL,SYMBOL,SYMBOL_SOURCE,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND" \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --symbol


# Perform VEP filtrations:
# Generate a grand VCF file with the manual filtration
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' ${data_path}/${project}/${subject}/${VEPOutputVcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence in missense_variant,frameshift_variant,inframe_deletion,coding_sequence_variant,stop_gained,inframe_insertion,splice_region_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is missense_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.missense_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.missense_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.missense_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is frameshift_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.frameshift_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.frameshift_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.frameshift_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is synonymous_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.synonymous_variant_only.query}


#######################################################################
# Variant loss against Parental
####################################################################### 

# grab difference as treatment variant gain from Parental
command_mem="100"
/home/luol2/lingqi_workspace/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${vcf_ctl} --discordance ${vcf_treatment} -O ${data_path}/${project}/${subject}/${sample_treatment}_vs_${sample_ctl}.vcf
# rename the vcf to Loss.xxxxxxx
tab=`echo -ne "\t"`;sed -i "s/\(CHROM.*FORMAT\).*/\1${tab}Loss.${sample_treatment}/" ${data_path}/${project}/${subject}/${sample_treatment}_vs_${sample_ctl}.vcf


module load singularity/3.1.1
VEPInputVcf=${sample_treatment}_vs_${sample_ctl}.vcf
VEPOutputVcf=${sample_treatment}_vs_${sample_ctl}.VEP.ann.vcf
VEPOutputStat=${sample_treatment}_vs_${sample_ctl}.VEP.ann.summary.html


#Perform VEP using local cache (Homo_Sapiens)
singularity run --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
vep -i /data/${VEPInputVcf} \
    -o /data/${VEPOutputVcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --stats_file /data/${VEPOutputStat}

#Perform VEP summarization after manual filtration
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' ${data_path}/${project}/${subject}/${VEPInputVcf} > ${data_path}/${project}/${subject}/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.vcf} 
singularity exec --bind ${data_path}/${project}/${subject}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.vcf} \
    -o /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.VEP.ann.vcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --stats_file /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.VEP.ann.summary.html}


# VEP output in tab format
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.vcf} \
    -o /data/${VEPInputVcf/.vcf/.AF0.05_BQ20_MQ50.VEP.ann.txt} \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Gene,CANONICAL,SYMBOL,SYMBOL_SOURCE,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND" \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --symbol


# Perform VEP filtrations:
# Generate a grand VCF file with the manual filtration
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' ${data_path}/${project}/${subject}/${VEPOutputVcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.query}


singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence in missense_variant,frameshift_variant,inframe_deletion,coding_sequence_variant,stop_gained,inframe_insertion,splice_region_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}//${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.non_synonymous_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is missense_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.missense_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.missense_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.missense_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is frameshift_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.frameshift_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.frameshift_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.frameshift_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}:/data \
/home/luol2/lingqi_workspace/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is synonymous_variant" | \
bcftools view -i 'INFO/DP > 10 & INFO/AF > 0.05 & INFO/BQ > 20 & INFO/MQ > 50' > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t%AF\t[%GT\t]\t\n" ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${VEPOutputVcf/.ann.vcf/.ann.AF0.05_BQ20_MQ50.synonymous_variant_only.query}
