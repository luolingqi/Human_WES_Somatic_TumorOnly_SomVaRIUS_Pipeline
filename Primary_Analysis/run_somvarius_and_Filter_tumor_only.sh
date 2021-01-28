#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 4 -R "rusage[mem=8]"
#BSUB -W 48:00

### Refer to  ~/lingqi_workspace/Projects/Mitesh/WGS/intel-gatk4-somatic-with-preprocessing/mutect2_nodocker.wdl

# prepare files that are required for SomVarIUS run

module load singularity/3.6.0
module load java/1.8.0_31

ref_fasta="Homo_sapiens_assembly38.fasta"
ref_dbsnps="dbSnp153Common_snvOnly.bed"
#/data/ldiaz/luol2/gatk_resource_bundle/hg38/dbSnp153Common.bed
tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"
command_mem="16000"

data_path=DATA_PATH
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

tumor_bam=${sample}.aligned.duplicates_marked.recalibrated.bam
tumor_vcf=${sample}.aligned.duplicates_marked.recalibrated.vcf

#tumor_bam=${sample}.aligned.duplicate_marked.sorted.bam.IndelRealigned.bam
#tumor_vcf=${sample}.aligned.duplicate_marked.sorted.bam.IndelRealigned.vcf

# index the bam
if [ ! -e ${data_path}/${project}/${subject}/${sample}/${tumor_bam/.bam/.bai} ];then
module load samtools/1.7
samtools index ${data_path}/${project}/${subject}/${sample}/${tumor_bam}
fi

# 1. Generate OXOG metrics:
/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" \
CollectSequencingArtifactMetrics \
-I ${data_path}/${project}/${subject}/${sample}/${tumor_bam} \
-O ${data_path}/${project}/${subject}/${sample}/${sample} \
--FILE_EXTENSION .txt \
-R /data/ldiaz/luol2/gatk_resource_bundle/hg38/${ref_fasta}

if [ 2 -eq 1 ]; then {
# SomVaRIUS calling
singularity run --bind ${data_path}/${project}/${subject}/${sample}:/data \
    --bind /data/ldiaz/luol2/gatk_resource_bundle/hg38:/ref \
    /home/luol2/lingqi_workspace/somvarius.sif SomVarIUS call_mutations \
    --bam /data/${tumor_bam} \
    --ref /ref/${ref_fasta} \
    --out /data/${tumor_vcf} \
    --germ_pos /ref/${ref_dbsnps/.bed/.pickle} \
    --dbsnp_bed /ref/${ref_dbsnps} \
    --dist /data/${tumor_bam/.bam/.betabinomial.dist.txt} \
    --min_pvalue 0.05 \
    --ref_filter
}
fi

## Rename the sample name in vcf
sed -i "s#/data/${sample}.aligned.duplicates_marked.recalibrated#${sample}#" ${data_path}/${project}/${subject}/${sample}/${tumor_vcf}

## Sort VCF with Picard
ref="/data/ldiaz/luol2/gatk_resource_bundle/hg38"

java -Xmx${command_mem}m -jar ${tool_path}/picard.jar \
SortVcf \
SEQUENCE_DICTIONARY=${ref}/Homo_sapiens_assembly38.dict \
OUTPUT=${data_path}/${project}/${subject}/${sample}/${tumor_vcf/.vcf/.targeted_sequencing.somvarius.tumor_only.sorted.vcf.gz} \
I=${data_path}/${project}/${subject}/${sample}/${tumor_vcf} \
CREATE_INDEX=true


#cp ${data_path}/${project}/${subject}/${sample}/${tumor_vcf}.stats ${data_path}/${project}/${subject}/${sample}/${tumor_vcf/.vcf/.targeted_sequencing.somvarius.tumor_only.sorted.vcf.gz}.stats


