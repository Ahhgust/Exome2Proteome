data="/home/augustw/src/Exome2Proteome/data"
bin="/home/augustw/src/Exome2Proteome/bin"
# above (data variable) needs to be re-defined at INSTALL 

# Files for BQSR (and VQSR):
indelFile="$data/SnpCalling/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
thousandGenomesSite="$data/SnpCalling/KG.phase1.snps.high_confidence.hg38.vcf.gz"  #original file: "1000G_phase3_v4_20130502.sites.vcf.gz", but a bug in GATK4 precludes filenames that start with integers. great.
dbSNP="$data/SnpCalling/dbsnp_138.hg38.vcf.gz"

# Additional VQSR file:
hapmapFile="$data/SnpCalling/hapmap_3.3.hg38.vcf.gz"

# everyone needs a human genome! This one has UCSC naming...
humanGenome="$data/HumanGenomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# a BED file for intervals targeted by Exome sequencing.
exomeBED="$data/ACE3_target.hg38.bed"
exomeBEDCommand="-L $exomeBED" # defined to empty string iff 3rd argument is "G"

# regex for BAMs
bamPattern=$1

if [ -z "$bamPattern" ]
then
    echo "$0 needs a *pattern* to work! E.g., $0 '*.bam', note the '', or foo.bam. note $0 foo.bam would generate only an individual's VCF"
    echo "Examples: \n$0 foo.bam 10 G (one person, 10 cores, whole genome (G) sequence)"
    echo "$0 'S*.bam' 20 G (all files that match S*.bam, 20 cores, whole genome sequence)"
    echo "$0 'S*.bam 20 (same as above, but whole EXOME sequence. Uses $data/ACE3_target.hg38.bed to define the exome target regions"
    echo "$0 'S*.bam (same as above, but uses the default core count (26). Uses $data/ACE3_target.hg38.bed to define the exome target regions"
    exit 1
fi

# The number of cores to use (with GNU's parallel)
# GATK will automatically add threads (so not helpful). So you may want to make this a small number
nCores=26

# used solely for picard's mark duplicates. good times that program. It automatically tries to use all cores. every time.
# this is set to a constant (2, which effectively oversubscribes the system by a factor of 2.... in practice, it works pretty well)
nCoresPicard=2

# ARG2 
if ! [ -z "$2" ]
then
    nCores=$2
    if ! [[ "$nCores" =~ ^[0-9]+$ ]]
    then
        echo "Optional argument $2 needs to be an integer (number of cores). E.g., $0 '*.bam' 10"
        exit 2
    fi
fi


# a BED file for intervals targeted by Exome sequencing.
exomeBED="$data/ACE3_target.hg38.bed"
exomeBEDCommand="-L $exomeBED"
if ! [ -z "$3" ]
then
    if [[ $3 -eq "G" ]]
    then
        exomeBEDCommand=""
    fi
fi


#TODO: Check for final VCFs. if found, exit!


parallelCommand="parallel --no-notice -j $nCores"
#parallelCommand="cat"
parallelCommandPicard="parallel --no-notice -j $nCoresPicard"


mkdir -p "bampreprocess"
# MARK DUPLICATES
for bam in $bamPattern; do
    if [ ! -f bampreprocess/`basename $bam .bam`.dedup.bam ]; then
        echo "PicardCommandLine MarkDuplicates QUIET=true I=$bam O=bampreprocess/`basename $bam .bam`.dedup.bam M=bampreprocess/`basename $bam .bam`.duplicateReadInfo"
    fi
done | $parallelCommandPicard || exit


# Begin BQSR: generate the list of co-variates
# this afaik can only be done using a single core, so it's a little slow (1 day)
# update: yes, in GATK4 this can be done one file at a time and thne GatherBQSRReports can be run.
# but no, I am unwilling to implement that....
if [ ! -f "bampreprocess/recal.csv" ]; then
    recalCommand="gatk BaseRecalibrator -R $humanGenome -O bampreprocess/recal.csv --known-sites $dbSNP --known-sites $indelFile --known-sites $thousandGenomesSite"
    for bam in $bamPattern; do
        dedupBam=bampreprocess/`basename $bam .bam`.dedup.bam
        recalCommand="$recalCommand -I $dedupBam"
    done
    # generate the table of covariates for BQSR
    $recalCommand || exit
fi


# Apply BQSR based on the table of covariates (recal.csv)
recalCommand="gatk ApplyBQSR -R $humanGenome --bqsr-recal-file bampreprocess/recal.csv"
for bam in $bamPattern; do
    dedupBam=bampreprocess/`basename $bam .bam`.dedup.bam
    outBam=bampreprocess/`basename $bam .bam`.dedup.bqsr.bam
    if [ ! -f $outBam ]; then
        echo "$recalCommand -I $dedupBam -O $outBam"
    fi
    
done | $parallelCommand || exit


mkdir -p "vcfout"

# GATK Haplotype caller
haplotypeCallerCommand="gatk HaplotypeCaller --dbsnp $dbSNP --genotyping-mode DISCOVERY  -R $humanGenome"
for bam in $bamPattern; do
    inBam=`basename $bam .bam`.dedup.bqsr.bam
    if [ ! -f bampreprocess/$inBam ]; then
        echo "Not good. Missing file: $inBam"
        exit
    fi
    haplotypeCallerCommand="$haplotypeCallerCommand -I bampreprocess/$inBam"
done

# apply multisample calling--
# Run in parallel across chromosomes
for chrom in {1..22} X Y M; do
    if [ ! -f vcfout/$vcfPrefix.$chrom.vcf.gz ]; then

        echo "$haplotypeCallerCommand -O vcfout/$vcfPrefix.$chrom.vcf.gz -L chr$chrom"
    fi
done | $parallelCommand || exit


# Modified from:
# https://software.broadinstitute.org/gatk/documentation/article?id=23216
# 
indelVQSRCommand="gatk --java-options \"-Xmx24g -Xms24g\" VariantRecalibrator --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -an QD -an SOR -mode INDEL --max-gaussians 4 -resource mills,known=false,training=true,truth=true,prior=12:$indelFile -resource dbsnp,known=true,training=false,truth=false,prior=2:$dbSNP -O vcfout/cohort_indels.recal --tranches-file vcfout/cohort_indels.tranches $exomeBEDCommand"

snpVQSRCommand="gatk --java-options \"-Xmx24g -Xms24g\" VariantRecalibrator --trust-all-polymorphic \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an InbreedingCoeff \
-mode SNP \
--max-gaussians 6 \
-resource hapmap,known=false,training=true,truth=true,prior=15:$hapmapFile \
-resource kg,known=false,training=true,truth=false,prior=10.:$thousandGenomesSite \
-resource dbsnp,known=true,training=false,truth=false,prior=7:$dbSNP \
-O vcfout/cohort_snps.recal \
--tranches-file vcfout/cohort_snps.tranches $exomeBEDCommand"

for chrom in {1..22} X Y M; do
    if [ -f "vcfout/$vcfPrefix.$chrom.vcf.gz" ]; then
        indelVQSRCommand="$indelVQSRCommand -V vcfout/$vcfPrefix.$chrom.vcf.gz"
        snpVQSRCommand="$snpVQSRCommand -V vcfout/$vcfPrefix.$chrom.vcf.gz"
    fi
done


if [ ! -f cohort_indels.recal.idx ]; then
    # generate the VQSR tables... (x2; once for SNPs, once for indels)
    echo -e "$snpVQSRCommand\n$indelVQSRCommand" | $parallelCommand || exit
fi

# apply the Gaussian mixture model
recalCommand="gatk --java-options \"-Xmx5g -Xms5g\" ApplyVQSR --recal-file vcfout/cohort_indels.recal \
--tranches-file vcfout/cohort_indels.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-mode INDEL"

# apply to INDELs
for chrom in {1..22} X Y M; do
    if [ ! -f $vcfPrefix.$chrom.INDEL.vqsr.vcf.gz ]; then
        echo "$recalCommand -V vcfout/$vcfPrefix.$chrom.vcf.gz -O vcfout/$vcfPrefix.$chrom.INDEL.vqsr.vcf.gz $exomeBEDCommand"
    fi
done | $parallelCommand || $exit

# and apply to SNPs
recalCommand="gatk --java-options \"-Xmx5g -Xms5g\" ApplyVQSR --recal-file vcfout/cohort_snps.recal \
--tranches-file vcfout/cohort_snps.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-mode SNP"

for chrom in {1..22} X Y M; do
    if [ ! -f $vcfPrefix.$chrom.recalibrated.vcf.gz ]; then
        echo "$recalCommand -O vcfout/$vcfPrefix.$chrom.recalibrated.vcf.gz -V vcfout/$vcfPrefix.$chrom.INDEL.vqsr.vcf.gz $exomeBEDCommand"
    fi
done | $parallelCommand || exit



