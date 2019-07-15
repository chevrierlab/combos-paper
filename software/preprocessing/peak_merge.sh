#Merges a bunch of .bed files into a unified peak set
#Takes a name for the comparison as the first input
#The remaining intputs are the names of the stimuli to be merged

#THIS PIECE OF CODE EXPECTS TO BE RUN IN THE 'OUT' DIRECTORY OF THE PROJECT IT'S WORKING ON

#setup input
name=$1
shift
module load bedtools
module load samtools

mkdir $name

#concatenate bed files
rm ${name}/${name}_concat.bed
for b in $@
do
	cat ${b}/${b}_peaks.narrowPeak >> ${name}/${name}_concat.bed
done

#sort concatenation
#samtools produces the order of the .bam, so the bed sorting matches the order of the .bam
#I'm not postive how helpful sorting is for merge, it seemed necessary
#for when I was using bedtools to count, but bedtools still needed too
#much memory
samtools view -H ${1}/${1}_sorted.bam | awk '$1 == "@SQ" {gsub ("SN:",""); gsub ("LN:", ""); print $2, $3}' | sed 's/ /\t/g' > ${name}/chr_order.txt
bedtools sort -faidx ${name}/chr_order.txt -i ${name}/${name}_concat.bed > ${name}/${name}_sorted_concat.bed

#merge intervals 
bedtools merge -sorted -i ${name}/${name}_sorted_concat.bed > ${name}/${name}_union.bed

#cleanup
rm ${name}/${name}_concat.bed ${name}/${name}_sorted_concat.bed ${name}/chr_order.txt
