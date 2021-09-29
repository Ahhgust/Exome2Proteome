# Assumptions: This assumes that the hs37d5 genome is used.
# This includes the chromosome names (1, MT, not chr1, chrM)
# as well as the 

# Get the requisite files:
# Bams
bamPattern=[SP]*.GRCh38.bam

# Write VCF files based on all bams
# the file prefix of these files is:
vcfPrefix=LLNL_HG38

# The number of cores to use (with GNU's parallel)
# GATK will automatically add threads (so not helpful). So you may want to make this a small number
nCores=26

# the FTP site to get the genomic data from
gatkFTP="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38"
# Files for BQSR (and VQSR):
indelFile="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
thousandGenomesSite="1000G_phase1.snps.high_confidence.hg38.vcf.gz"  #"1000G_phase3_v4_20130502.sites.vcf.gz"
dbSNP="dbsnp_138.hg38.vcf.gz"

# Additional VQSR file:
hapmapFile="hapmap_3.3.hg38.vcf.gz"

parallelCommand="parallel --no-notice -j $nCores"
parallelCommandJava="parallel --no-notice -j 6"

### Local resources
# Refreence genome
humanGenome="/home/data/LabData/HumanGenomes/GRCh38_Recommended/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Same genome using ENSEMBL identifiers...
humanGenomeEnsembl="/home/data/LabData/HumanGenomes/GRCh38_Recommended/GRCh38_full_analysis_set_plus_decoy_hla.ensembl.fa"

# Exome BED:
exomeBED="ACE3_target.hg38.bed"

# Protein coding genes (fasta) from Ensembl
ensemblPep="/home/data/LabData/HumanGenomes/Ensembl/Ensembl_101/Homo_sapiens.GRCh38.pep.all.fa.gz"

# boolean; do we need to convert from UCSC (chr1) to Ensembl (1) chromosome naming?
convertUCSC2ensembl=true


# The Broad FTP site has some strict limits on the number of concurrent users
# This tries to download a file from this site (or any other site)
# and if the download fails, it retries a set number of times 
function downloadHelper {

    # only download once...
    if [ -f $2 ]; then
        return 2
    fi
    if [ -f `basename $2 .gz` ]; then
        return 3
    fi

    sleepLength=1
    echo "Attempting to download $2"
    for i in {1..20}; do
        echo "Attempt $i"
        wget --quiet "$1/$2"
        if [ $? -eq 0 ]; then
            return 1
        fi
        sleep $sleepLength
        sleepLength=$((sleepLength+1))
    done

    echo "Too many download attempts for file: $1/$2"
    exit 
}

# download files and indexes.
downloadHelper $gatkFTP $dbSNP
downloadHelper $gatkFTP `basename $dbSNP`.tbi

downloadHelper $gatkFTP $indelFile
downloadHelper $gatkFTP `basename $indelFile`.tbi


downloadHelper $gatkFTP $thousandGenomesSite
downloadHelper $gatkFTP `basename $thousandGenomesSite`.tbi

downloadHelper $gatkFTP $hapmapFile
downloadHelper $gatkFTP `basename $hapmapFile`.tbi



# remove the .gz file extensions (downloadHelper also gunzips..)
#dbSNP=`basename $dbSNP .gz`
#indelFile=`basename $indelFile .gz`
#thousandGenomesSite=`basename $thousandGenomesSite .gz`
#hapmapFile=`basename $hapmapFile .gz`

## Joy!
# There's a bug in GATK; it can't recognize file names that begin with a number...
# (when specified as a resource)
# hence the 1000 genomes file gets a new name
if [ ! -f "KG.phase1.snps.high_confidence.hg38.vcf.gz" ]; then
    ln -s $thousandGenomesSite KG.phase1.snps.high_confidence.hg38.vcf.gz
    ln -s $thousandGenomesSite.gz.tbi KG.phase1.snps.high_confidence.hg38.vcf.gz.tbi
fi

thousandGenomesSite="KG.phase1.snps.high_confidence.hg38.vcf.gz"

# check to see if the appropriate programs are found...
# picard
PicardCommandLine -h 2>/dev/null
if [ $? -ne 1 ]; then # don't ask me why it returns 1.
    echo "Please install Picard!"
    exit
fi

# GNU's parallel
parallel -h >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install GNU parallel!"
    exit
fi

# GATK (4.*)
gatk --help >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install GATK!"
    exit
fi

#what's hap! (physical phasing)
whatshap --help >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install whatshap!!"
    exit
fi


bcftools -h >/dev/null
if [ $? -ne 0 ]; then # check for the correct return value (0)
    echo "Please install bcftools!"
    exit
fi

tabix -h 2>/dev/null
if [ $? -ne 1 ]; then # check for the correct return value (1 in this case... 'cause they're annoying)
    echo "Please install tabix!!"
    exit
fi

shapeit4 &> /dev/null
if [ $? -ne 1 ]; then # check for the return val
    echo "Please install shapeit4!!"
    exit
fi

# No auto-download for these files... they're too big

if [ ! -d 1000Genomes ]; then
    echo "Please download the 1000 genomes data and place it in a directory 1000Genomes. I'm looking for the 20181129  genotypes vcfs..."
    echo "One such site is: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    exit
fi


# Mark PCR duplicates
# in parallel
for bam in $bamPattern; do
    if [ ! -f `basename $bam .bam`.dedup.bam ]; then
        echo "PicardCommandLine MarkDuplicates QUIET=true I=$bam O=`basename $bam .bam`.dedup.bam M=`basename $bam .bam`.duplicateReadInfo"
    fi
done | $parallelCommandJava || exit


# Begin BQSR: generate the list of co-variates
# this afaik can only be done using a single core, so it's a little slow (1 day)
# update: yes, in GATK4 this can be done one file at a time and thne GatherBQSRReports can be run.
if [ ! -f recal.csv ]; then
    recalCommand="gatk BaseRecalibrator -R $humanGenome -O recal.csv --known-sites $dbSNP --known-sites $indelFile --known-sites $thousandGenomesSite"
    for bam in $bamPattern; do
        dedupBam=`basename $bam .bam`.dedup.bam
        recalCommand="$recalCommand -I $dedupBam"
    done
    # generate the table of covariates for BQSR
    $recalCommand || exit
fi


# Apply BQSR based on the table of covariates (recal.csv)
recalCommand="gatk ApplyBQSR -R $humanGenome --bqsr-recal-file recal.csv"
for bam in $bamPattern; do
    dedupBam=`basename $bam .bam`.dedup.bam
    outBam=`basename $bam .bam`.dedup.bqsr.bam
    if [ ! -f $outBam ]; then
        echo "$recalCommand -I $dedupBam -O $outBam"
    fi
    
done | $parallelCommand || exit

# GATK Haplotype caller
haplotypeCallerCommand="gatk HaplotypeCaller --dbsnp $dbSNP --genotyping-mode DISCOVERY  -R $humanGenome"
for bam in $bamPattern; do
    inBam=`basename $bam .bam`.dedup.bqsr.bam
    if [ ! -f $inBam ]; then
        echo "Not good. Missing file: $inBam"
        exit
    fi
    haplotypeCallerCommand="$haplotypeCallerCommand -I $inBam"
done

# apply multisample calling--
# Run in parallel across chromosomes
for chrom in {1..22} X Y M; do
    if [ ! -f $vcfPrefix.$chrom.vcf.gz ]; then

        echo "$haplotypeCallerCommand -O $vcfPrefix.$chrom.vcf.gz -L chr$chrom"
    fi
done | $parallelCommand || exit


# Modified from:
# https://software.broadinstitute.org/gatk/documentation/article?id=23216
# 
indelVQSRCommand="gatk --java-options \"-Xmx24g -Xms24g\" VariantRecalibrator --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -an QD -an SOR -mode INDEL --max-gaussians 4 -resource mills,known=false,training=true,truth=true,prior=12:$indelFile -resource dbsnp,known=true,training=false,truth=false,prior=2:$dbSNP -O cohort_indels.recal --tranches-file cohort_indels.tranches -L $exomeBED"

snpVQSRCommand="gatk --java-options \"-Xmx24g -Xms24g\" VariantRecalibrator --trust-all-polymorphic \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an InbreedingCoeff \
-mode SNP \
--max-gaussians 6 \
-resource hapmap,known=false,training=true,truth=true,prior=15:$hapmapFile \
-resource kg,known=false,training=true,truth=false,prior=10.:$thousandGenomesSite \
-resource dbsnp,known=true,training=false,truth=false,prior=7:$dbSNP \
-O cohort_snps.recal \
--tranches-file cohort_snps.tranches \
-L $exomeBED"

for chrom in {1..22} X Y M; do
    if [ -f "$vcfPrefix.$chrom.vcf.gz" ]; then
        indelVQSRCommand="$indelVQSRCommand -V $vcfPrefix.$chrom.vcf.gz"
        snpVQSRCommand="$snpVQSRCommand -V $vcfPrefix.$chrom.vcf.gz"
    fi
done


if [ ! -f cohort_indels.recal.idx ]; then
    # generate the VQSR tables... (x2; once for SNPs, once for indels)
    echo -e "$snpVQSRCommand\n$indelVQSRCommand" | $parallelCommand || exit
fi

# apply the Gaussian mixture model
recalCommand="gatk --java-options \"-Xmx5g -Xms5g\" ApplyVQSR --recal-file cohort_indels.recal \
--tranches-file cohort_indels.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-mode INDEL"

# apply to INDELs
for chrom in {1..22} X Y M; do
    if [ ! -f $vcfPrefix.$chrom.INDEL.vqsr.vcf.gz ]; then
        echo "$recalCommand -V $vcfPrefix.$chrom.vcf.gz -O $vcfPrefix.$chrom.INDEL.vqsr.vcf.gz -L $exomeBED"
    fi
done | $parallelCommand || $exit

# and apply to SNPs
recalCommand="gatk --java-options \"-Xmx5g -Xms5g\" ApplyVQSR --recal-file cohort_snps.recal \
--tranches-file cohort_snps.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-mode SNP"

for chrom in {1..22} X Y M; do
    if [ ! -f $vcfPrefix.$chrom.recalibrated.vcf.gz ]; then
        echo "$recalCommand -O $vcfPrefix.$chrom.recalibrated.vcf.gz -V $vcfPrefix.$chrom.INDEL.vqsr.vcf.gz -L $exomeBED"
    fi
done | $parallelCommand || exit

### PHASING BEGINS HERE

# start with physical pre-phasing

# phases indels; samples 20 reads to determine phasing
# scales as ~2^20, so 20 is a practical upper-limit
physicalPhaseCommand="whatshap phase --indels --max-coverage 20"
# note no Y and MT
# though I have a morbid curiousity as to what it would do.
# at least the Y could have some value...
for chrom in {1..22} X; do
    if [ ! -f $vcfPrefix.$chrom.recalibrated.whatshap.vcf.gz ]; then
        echo "$physicalPhaseCommand -o $vcfPrefix.$chrom.recalibrated.whatshap.vcf.gz $vcfPrefix.$chrom.recalibrated.vcf.gz *.dedup.bqsr.bam"
    fi
done | $parallelCommand || exit

# and then do statistical phasing

# first, we need to convert the vcf.gz files into index bcf files.
# BCF
for chrom in {1..22} X; do
    if [ ! -f $vcfPrefix.$chrom.recalibrated.whatshap.bcf ]; then
        echo "bcftools concat -O b -o $vcfPrefix.$chrom.recalibrated.whatshap.bcf $vcfPrefix.$chrom.recalibrated.whatshap.vcf.gz "
    fi
done | $parallelCommand || exit

for chrom in {1..22} X; do
    if [ ! -f $vcfPrefix.$chrom.recalibrated.whatshap.bcf.csi ]; then
        echo "bcftools index $vcfPrefix.$chrom.recalibrated.whatshap.bcf "
    fi
done | $parallelCommand || exit


# and generate the sample-specific fastas...
ids=`zcat $vcfPrefix.22.recalibrated.whatshap.vcf.gz | fgrep '#CHROM' | cut -f10- | tr '\t' '\n'`


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
        if [ ! -f $vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf.csi ]; then
            # was with -R exomeBED; not necessary. may be bug with multisample...
# wrong. 1KG genotype data are in UCSC format. need to punt.
            #            if [ convertUCSC2ensembl ]; then
 #               echo "bcftools annotate --rename-chrs chr_name_conv.txt -O u $vcfPrefix.$chrom.recalibrated.whatshap.bcf | bcftools view -l 9 -f PASS -s $id -I -a -O b -o $vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf /dev/stdin && bcftools index $vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf"
  #          else
                echo "bcftools view -l 9 -f PASS -s $id -I -a -O b -o $vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf $vcfPrefix.$chrom.recalibrated.whatshap.bcf && bcftools index $vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf"
   #         fi
        fi
    done
done | $parallelCommand || exit


for chrom in {1..22} X; do
    # and phase!
    # use a reference panel and physical phasing
    # sites (matching on position+alt allele) not in the reference panel are dropped.
    for id in $ids ; do
        if [ ! -f $vcfPrefix.$chrom.$id.recalibrated.whatshap.shapeit4.bcf.csi ]; then

            seed=$chrom
            if [ $seed == 'X' ]; then
                seed='23'
            fi
           echo "shapeit4 --input $vcfPrefix.$chrom.$id.recalibrated.whatshap.bcf --map ./GeneticMaps/chr$chrom.b38.gmap.ucsc.gz --region chr$chrom --use-PS 0.0001 --log $vcfPrefix.$chrom.$id.shapeit4Log --sequencing --seed $seed --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --output $vcfPrefix.$chrom.$id.recalibrated.whatshap.shapeit4.bcf -H 1000Genomes/ALL.chr$chrom.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz && bcftools index $vcfPrefix.$chrom.$id.recalibrated.whatshap.shapeit4.bcf "
        fi
    done
done  | $parallelCommand || exit

for chrom in {1..22} X; do
    # and merge the single-sample phased BCFs (with some non-panel alleles dropped)
    # and harmonize; add in the dropped genotypes.
    # -f: the shapeit4 bcfs have the filter=='.'
    if [ ! -f $vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz ]; then
        # was with -R exomeBED; not necessary. and bug with multisample...
        echo "bcftools merge -f PASS,. -0 -O v --force-samples -m all $vcfPrefix.$chrom.?*recalibrated.whatshap.shapeit4.bcf $vcfPrefix.$chrom.recalibrated.whatshap.bcf | python3 harmonizePhasingMerge.py $vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz"
    fi
    
done | $parallelCommand || exit


#TODO: Add step to convert from UCSC to Ensembl nomenclature here
# ALSO UPDATE GFFs

for chrom in {1..22}; do
    # run bcftools csq
    # the grep '>.*|' serves to just keep the variants that change amino acids...
    if [ ! -f $vcfPrefix.$chrom.ProteinCodingConsequences.txt.gz ]; then
        if [ convertUCSC2ensembl ]; then
            echo "bcftools annotate --rename-chrs chr_name_conv.txt -O u $vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz | bcftools csq -f $humanGenomeEnsembl -g Ensembl_101/Homo_sapiens.GRCh38.101.chromosome.$chrom.gff3.gz -O t -p a /dev/stdin |  grep '>.*|' | gzip -9 > $vcfPrefix.$chrom.ProteinCodingConsequences.txt.gz" 
        else
            echo "bcftools csq -f $humanGenome -g Ensembl_101/Homo_sapiens.GRCh38.101.chromosome.$chrom.gff3.gz -O t -p a $vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz |  grep '>.*|' | gzip -9 > $vcfPrefix.$chrom.ProteinCodingConsequences.txt.gz"
        fi
    fi
done | $parallelCommand || exit


for id in $ids; do
    if [ ! -f $id.fa.gz ]; then
        echo "./gvp2null.py -k -j $id  -c '$vcfPrefix*ProteinCodingConsequences.txt.gz' -p $ensemblPep | gzip -9 > $id.fa.gz"
    fi
done | $parallelCommand || exit

# and the same data for all individuals combined...
if [ ! -f $vcfPrefix.fa.gz ]; then
    echo "./gvp2null.py -d $vcfPrefix.LUT.gz -c '$vcfPrefix*ProteinCodingConsequences.txt.gz' -p $ensemblPep | gzip -9 > $vcfPrefix.fa.gz"
fi | $parallelCommand || exit



