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
**Each BAM must also be aligned to the GRCh38 reference genome. This means UCSC-style chromosome naming (e.g., chr1, not 1). Only the canonical autosomes (chr1-chr22 are used)**

For whole exome data:
```
  bin/bam2vcf.sh '*.bam' 23 &> variantcalling.outerr 
```
Uses up to 23 cores*
<br>
<br>
(actually that's a lie, Picard auto-magically uses way more than that, but it's not a parameter you can set)

For whole genome data:
```
  bin/bam2vcf.sh '*.bam' 26 G &> variantcalling.outerr 
```
Which (you guessed it) uses up to 26 cores.

Unphased VCFs should be in the directory ```vcfout/```
And have names like ```E2P.1.recalibrated.vcf.gz```
(meaning exome 2 proteome, chromosome 1, VQSR recalibrated)

## Unphased VCFs -> phased VCFs
Does physical phasing (whatshap) followed by statistical phasing (shapeit4). Recovers sites lost (dropped by shapeit4) in the statistical phasing step.
Type:
```
  bin/vcf2phase.sh 23 &> phasing.outerr 
```
This uses up to 23 cores. The final vcfs are written to ``finalvcf/``` (one vcf per chromosome)
A directory called ```phasedvcf/``` is also made; it has temporary files (each chromosome x each person) and is only necessary for debugging.

## Phased VCF -> Proteomes in fasta

```
  bin/phasedvcf2prot.sh 23 &> proteome.outerr 
```


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



