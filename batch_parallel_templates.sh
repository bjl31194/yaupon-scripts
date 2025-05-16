## make list of names of individuals
INDS=($(for i in /home/data/wgs_raw/*R1.fastq.gz; do echo $(basename ${i%.R*}); done))

## run a script on list of individuals iteratively
for IND in ${INDS[@]};
do
	# declare variables
	FORWARD=/home/data/wgs_raw/${IND}.R1.fastq.gz
	REVERSE=/home/data/wgs_raw/${IND}.R2.fastq.gz
	OUTPUT=~/align/${IND}_sort.bam

    # then align and sort
	echo "Aligning $IND with bwa"
	bwa-mem2 mem -t 4 $REF $FORWARD \
	$REVERSE | samtools view -b | \
	samtools sort -T ${IND} > $OUTPUT

done

## make a file with individual names 
for i in /home/data/wgs_raw/*R1.fastq.gz; do echo $(basename ${i%.R*}); done > inds

## script template to run with parallel called "parallel_align.sh"
# align a single individual
REF=~/reference/P_nyererei_v2.fasta

# declare variables
IND=$1
FORWARD=/home/data/wgs_raw/${IND}.R1.fastq.gz
REVERSE=/home/data/wgs_raw/${IND}.R2.fastq.gz
OUTPUT=~/align/${IND}_sort.bam

# then align and sort
echo "Aligning $IND with bwa"
bwa-mem2 mem -t 4 $REF $FORWARD \
$REVERSE | samtools view -b | \
samtools sort -T ${IND} > $OUTPUT

## run script "parallel_align.sh" in parallel using "inds" file
parallel `sh parallel_align.sh {}` :::: inds