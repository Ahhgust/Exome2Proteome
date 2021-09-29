
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


if ! [ -d "finalvcf" ]
then
    echo "Expected directory: finalvcf/ is not found. Are you calling this program from the right directory?"
    echo "Usage: (run bam2vcf.sh first!)\n$0 FOO 30"
    echo "Which will generate multi-sample phased VCFs with named phasedvcf/FOO.*bcf, using up to 30 cores at a time"
    exit 1
fi

humanGenomeEnsembl="$data/HumanGenomes/GRCh38_full_analysis_set_plus_decoy_hla.ensembl.fa"
parallelCommand="parallel --no-notice -j $nCores"

mkdir -p "proteome"

for chrom in {1..22}; do
    # run bcftools csq
    # the grep '>.*|' serves to just keep the variants that change amino acids...
    if [ ! -f proteome/$vcfPrefix.$chrom.ProteinCodingConsequences.txt.gz ]; then
        echo "bcftools annotate --rename-chrs $data/chr_name_conv.txt -O u finalvcf/$vcfPrefix.$chrom.recalibrated.whatshap.shapeit4.vcf.gz | bcftools csq -f $humanGenomeEnsembl -g $data/PredictingProteins/Ensembl_101/Homo_sapiens.GRCh38.101.chromosome.$chrom.gff3.gz -O t -p a /dev/stdin |  grep '>.*|' | gzip -9 > proteome/$vcfPrefix.$chrom.ProteinCodingConsequences.txt.gz" 
    fi
done | $parallelCommand || exit


for id in $ids; do
    if [ ! -f proteome/$id.fa.gz ]; then
        echo "$bin/gvp2null.py -k -j $id  -c 'proteome/$vcfPrefix*ProteinCodingConsequences.txt.gz' -p $data/PredictingProteins/Ensembl_101/Homo_sapiens.GRCh38.pep.all.fa.gz | gzip -9 > proteome/$id.fa.gz"
    fi
done | $parallelCommand || exit

# and the same data for all individuals combined...
if [ ! -f proteome/$vcfPrefix.fa.gz ]; then
    echo "$bin/gvp2null.py -d proteome/$vcfPrefix.LUT.gz -c 'proteome/$vcfPrefix*ProteinCodingConsequences.txt.gz' -p $data/PredictingProteins/Ensembl_101/Homo_sapiens.GRCh38.pep.all.fa.gz | gzip -9 > proteome/$vcfPrefix.fa.gz"
fi | $parallelCommand || exit




