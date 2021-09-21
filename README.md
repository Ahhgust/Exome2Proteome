#

How to convert whole genomes (or exomes) into proteomes!
The general flow of information is to take raw bams (Mapped to the GRCh38 reference genome) and
* Call variants (GATK4 currently. There is a boatload of steps involved here. We plan on moving to DeepVariant)
  * process bams by
    * marking duplicates
    * recalibration base qualities (BQSR)
  * Call snps (GATK's haplotypecaller)
    * Identify true variants
      * GATK's VQSR
* Phase variants
  * Physically (WhatsHap)
    * This is important for SNP variants that are close to each other!
  * Statistically (shapeit4)
* Predict *haploid* protein sequences (bcftools csq)
  * Produce fasta sequences from these predictions (custom)
  * An index these fasta files with a suffix array (protengine)


# Installation (*nix only)
This involves downloading a boatload of genomic data.

```
git clone https://github.com/Ahhgust/Exome2Proteome.git && cd Exome2Proteome

wget --quiet -O data.tbz 'https://www.dropbox.com/s/jt9doz7ixmovhik/DataBundle_Genomes2Proteome.tbz?dl=1'

tar -xf data.tbz && rm data.tbz

(echo data=\"$PWD/data\"; echo bin=\"$PWD/bin\" ; tail -n +3 bin/bam2vcf.sh) > foo && mv -f foo bin/bam2vcf.sh
(echo data=\"$PWD/data\"; echo bin=\"$PWD/bin\" ; tail -n +3 bin/vcf2prot.sh) > foo && mv -f foo bin/vcf2prot.sh

```

# Quick start

# Depedency check






# Data input

For bam2vcf:

I assume that BAM files have been created; *one BAM per individual.*
BAM files should be mapped but unprocessed.
These files need to be aligned to the GRCh38 reference genome. In particular,
they need to have chromosome names like chr22 (not 22!).
The scripts only consider the canonical autosomes: (chr1:chr22). ChrX is a work in progress.
<br>
The output from bam2vcf are VCF phased VCF files, with those variants that were deemed "unphasable"
retained.


# Dependencies (programs) (way too many! *nix only)

All dependencies must be in your $PATH (ie, type any of the following , and it works)

* gatk              (v4. Must have HaplotypeCaller))
* PicardCommandLine (picard)
* bcftools
* shapeit4
* whatshap
* bwa
* samtools

# Dependencies (data)
There are so many data dependencies, and they are all incredibly particular.
As such, I bundled them for you. Apologies for any redundancies...






