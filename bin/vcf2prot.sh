data="/home/augustw/src/Exome2Proteome/data"
bin="/home/augustw/src/Exome2Proteome/bin"
# above (data variable) needs to be re-defined at INSTALL 

# The number of cores to use (with GNU's parallel)
# GATK will automatically add threads (so not helpful). So you may want to make this a small number
nCores=26

vcfPrefix="HG38"
if ! [ -z "$1" ]
then
    vcfPrefix=$1
fi


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


if ! [ -d "bampreprocess" ]
then
    echo "Expected directory: bampreprocess/ is not found. Are you calling this program from the right directory?"
    echo "Usage: (run bam2vcf.sh first!)\n$0 FOO 30"
    echo "Which will generate multi-sample phased VCFs with named phasedvcf/FOO.*bcf, using up to 30 cores at a time"
    exit 1
fi

if ! [ -d  "vcfout" ]
then
    echo "Expected directory: vcfout/ is not found. Are you calling this program from the right directory?"
    echo "Usage: (run bam2vcf.sh first!)\n$0 FOO 30"
    echo "Which will generate multi-sample phased VCFs with named phasedvcf/FOO.*bcf, using up to 30 cores at a time"
    exit 1
fi

parallelCommand="parallel --no-notice -j $nCores"

mkdir -p "phasedvcf"
# phases indels; samples 20 reads to determine phasing
# scales as ~2^20, so 20 is a practical upper-limit
physicalPhaseCommand="whatshap phase --indels --max-coverage 20"
# note no Y and MT
# though I have a morbid curiousity as to what it would do.
# at least the Y could have some value...
# modified: now creates BCF directly
for chrom in {1..22} X; do
    if [ ! -f phasedvcf/$vcfPrefix.$chrom.recalibrated.whatshap.bcf ]; then
        echo "$physicalPhaseCommand vcfout/$vcfPrefix.$chrom.recalibrated.vcf.gz bampreprocess/*.dedup.bqsr.bam | bcftools view -l 9 --threads 3 -O b -o phasedvcf/$vcfPrefix.$chrom.recalibrated.whatshap.bcf /dev/stdin"
    fi
done | $parallelCommand || exit

for chrom in {1..22} X; do
    if [ ! -f phasedvcf/$vcfPrefix.$chrom.recalibrated.whatshap.bcf.csi ]; then
        echo "bcftools index phasedvcf/$vcfPrefix.$chrom.recalibrated.whatshap.bcf"
    fi
done | $parallelCommand || exit




