data="/home/augustw/src/Exome2Proteome/data"
bin="/home/augustw/src/Exome2Proteome/bin"
# above (data variable) needs to be re-defined at INSTALL 

# The number of cores to use (with GNU's parallel)
# GATK will automatically add threads (so not helpful). So you may want to make this a small number
nCores=26

vcfPrefix="E2P"
# ARG1 
if ! [ -z "$1" ]
then
    nCores=$1
    if ! [[ "$nCores" =~ ^[0-9]+$ ]]
    then
        echo "Optional argument $1 needs to be an integer (number of cores). E.g., $0 10"
        exit 1
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

#ids=`zcat $vcfPrefix.22.recalibrated.whatshap.vcf.gz | fgrep '#CHROM' | cut -f10- | tr '\t' '\n'`
ids=`bcftools query -l phasedvcf/$vcfPrefix.22.recalibrated.whatshap.bcf`


#TODO:
# split bcf into 1 bcf per person
# phase that with the 1KG as a reference
# and merge/harmonize into multisample callign
# need to specify the phased files first
# bcftools merge -f PASS -0  --force-samples -m all  16.test.PR17.bcf 16.test.Pr1.bcf PR17.whatshap.bcf Pr1.whatshap.bcf -O v  | python3 harmonizePhasingMerge.py  harmony.vcf.gz

# restrict sites to just the exome; merge the 1000 Genomes data with the data we have
# since it's all exome sequencnig, call the non-called sites as 0/0 (imperfect but adequate assumption)
# then split out multiallelic sites into multiple rows of biallelic sites (biallelics only for phasing)
# and then sort (shouldn've have to do this, but you do)
for chrom in {1..22} X; do
    # make a sample-specific BCF (and chromosome specific)
    # restrict to just the exome
    # and VQSR PASS sites
    # and extra alleles get dropped (if not present in person)
    # and if this process executes w/o error, make the index as well
    for id in $ids ; do
        if [ ! -f phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf.csi ]; then
            echo "bcftools view -l 9 -f PASS -s $id -I -a -O b -o phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf phasedvcf/$vcfPrefix.$chrom.recalibrated.whatshap.bcf && bcftools index phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf"
        fi
    done
done | $parallelCommand || exit


# do statistical phasing on each individual, each chromosome (in parallel)
for chrom in {1..22} X; do
    # and phase!
    # use a reference panel and physical phasing
    # sites (matching on position+alt allele) not in the reference panel are dropped.
    for id in $ids ; do
        if [ ! -f phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.shapeit4.bcf.csi ]; then
            seed=$chrom
            if [ $seed == 'X' ]; then
                seed='23'
            fi
           echo "shapeit4 --input phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf --map $data/Phasing/GeneticMaps/chr$chrom.b38.gmap.ucsc.gz --region chr$chrom --use-PS 0.0001 --log phasedvcf/$vcfPrefix.$chrom.$id.shapeit4Log --sequencing --seed $seed --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --output phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.shapeit4.bcf -H $data/Phasing/ALL.chr$chrom.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz && bcftools index phasedvcf/$vcfPrefix.$chrom.$id.recalibrated.whatshap.shapeit4.bcf "
        fi
        
    done
done  | $parallelCommand || exit

mkdir -p "finalvcf"

for chrom in {1..22} X; do
    # and merge the single-sample phased BCFs (with some non-panel alleles dropped)
    # and harmonize; add in the dropped genotypes.
    # -f: the shapeit4 bcfs have the filter=='.'
    if [ ! -f "finalvcf/$vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz" ]; then
        # was with -R exomeBED; not necessary. and bug with multisample...
        echo "bcftools merge -f PASS,. -0 -O v --force-samples -m all phasedvcf/$vcfPrefix.$chrom.?*recalibrated.whatshap.shapeit4.bcf phasedvcf/$vcfPrefix.$chrom.recalibrated.whatshap.bcf | python3 $bin/harmonizePhasingMerge.py finalvcf/$vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz"
    fi
    
done | $parallelCommand || exit


