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
This involves downloading a boatload of genomic data (wget).
As well as a poor-man's Make to embed where the data live.
```
git clone https://github.com/Ahhgust/Exome2Proteome.git && cd Exome2Proteome
wget --quiet -O data.tbz 'https://www.dropbox.com/s/jt9doz7ixmovhik/DataBundle_Genomes2Proteome.tbz?dl=1'
tar -xf data.tbz && rm data.tbz
(echo data=\"$PWD/data\"; echo bin=\"$PWD/bin\" ; tail -n +3 bin/bam2vcf.sh) > foo && mv -f foo bin/bam2vcf.sh
(echo data=\"$PWD/data\"; echo bin=\"$PWD/bin\" ; tail -n +3 bin/vcf2prot.sh) > foo && mv -f foo bin/vcf2prot.sh
chmod +x bin/*
```

# Quick start

## Dependency check
Run:
```
  bin/dependency_check.sh
```
This is a case of silence is golden-- if the script prints nothing and exits cleanly, then all of the *programs* needed to convert exomes to protomes have been installed.
I leave it up to the user to ensure that the programs are reasonably up-to-date (hint hint you may need to update bcftools)

## Bams -> unphased VCFs
Convert bam to (unphased) VCF. GATK4-style.
Assumptions: you're in a directory with 1+ bams in it (and nothing else).
**Each BAM file must correspond to a single sample. Read-groups must be set and valid**

For whole exome data:
```
  bin/bam2vcf.sh '*.bam' 23 &> variantcalling.outerr &
```
Uses up to 23 cores*
<br>
<br>
(actually that's a lie, Picard auto-magically uses way more than that, but it's not a parameter you can set)

For whole genome data:
```
  bin/bam2vcf.sh '*.bam' 26 G &> variantcalling.outerr &
```

Unphased VCFs should be in the directory ```vcfout/```

## Unphased VCFs -> phased VCFs
Does physical phasing followed by statistical phasing. Recovers sites lost (dropped by shapeit4) in the statistical phasing step.

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
* python3
* R

Additionally, the following R packages are needed for generating an in silico trypsin digest:
* tools
* R.utils
* optparse
* cleaver
* Peptides
* stringr
* stringi
* seqinr
* yaml
* data.table



