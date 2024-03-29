
# define path for reference genome and index
bwaIndex=/path/to/bwa_index
gff3=/path/to /annotation.gff3
# Create directory for filtered read and count matrix
mkdir -p filtered_read
mkdir -p htseq_count

# Looping over read1
for R1 in *_R1_*.gz
do
    # replacing _R1_ with _R2_ for reverse read. please look at your pattern
    R2=(${R1//_R1_/_R2_})
     # you can print with echo command to see whether it is working or not.
    #echo $R1 $R2 

    #create a variable B for base name for reading convenience
   # which only captures the sample information
    B=${R1::-16}
   # Here B trims all the last 16 character from read file. You can choose your own number
    fastp  -i $R1 -I $R2  --detect_adapter_for_pe  -o filtered_read/$R1 -O filtered_read/$R2  -h filtered_read/$B.html
    bwa mem -t 8 $bwaIndex filtered_read/$R1 filtered_read/$R2 > ${B}.sam
    samtools view -Sb ${B}.sam > ${B}.unsorted.bam
    samtools sort  -o ${B}_sorted.bam ${B}.unsorted.bam
    samtools index ${B}.bam
    htseq-count --format=bam --stranded=no --type=gene --order=pos --idattr=ID ${B}_sorted.bam  $gff3  > htseq_count/${B}.txt

done
