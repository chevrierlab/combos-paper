#!/bin/bash
#SBATCH --job-name=atacseq_pipeline
#SBATCH --output=/project2/nchevrier/agtree/combos/atacseq/logs/pipeline.out
#SBATCH --error=/project2/nchevrier/agtree/combos/atacseq/logs/pipeline.err
#SBATCH --tasks-per-node 8  #Number of cores 
#SBATCH --time 24:00:00 
#SBATCH --partition broadwl  #Partition to submit to (other option: general)
#SBATCH --mem=16384 #Memory per node in MB (see also --mem-per-cpu) 


#read in arguments. Data folder first, then name, then flags (mostly for debugging)
#available flags:
	#-nc: keeps the pipeline from cleaning up temporary files
	#-skip [skip] : skips the pipeline to a midpoint
dir=$1
shift
name=$1
shift

#this is some code I cribbed from an example for how to make flags in bash
#each iteration of the loop looks for a flag, and then shifts the inputs one
while test $# -gt 0; do
    case "$1" in 
        -nc)
            nc=nc
            shift
            ;;
        -skip)
            shift
            if test $# -gt 0; then
                export skip=$1
            else
                echo "no skip point specified"
                exit 1
            fi
            shift
            ;;
        *)
            break
            ;;
    esac
done

#housekeeping
olddir=`pwd`
cd $dir
mkdir out/${name}

#this skip methodology is used to allow someone to skip over the first portion of a pipeline
#(in case it failed at a point, and they don't want to repeat analysis)
case "$skip" in 
    "")
        #trim
        ;&
    bowtie)
        #bowtie aligns reads. I decided to use a hard link to the index -- something to potentially reevaluate
        module load bowtie2
        bowtie2 -X 2000 -p 8 -x /project2/nchevrier/agtree/indices/mm10_bowtie_index/mm10 -1 fastq/${name}_R1.fastq.gz -2 fastq/${name}_R2.fastq.gz -S out/${name}/${name}.sam 
        ;&
    samtools)
        #samtools
        module load samtools
        samtools view -bS out/${name}/${name}.sam | samtools sort -o out/${name}/${name}_sorted.bam
        if [[ -z "$nc" ]] 
        then
            rm out/${name}/${name}.sam 
        fi
        ;&
    igv)
        #igv
        module load igv
        igvtools count -w 5 out/${name}/${name}_sorted.bam out/${name}/${name}.tdf mm10
        ;&
    macs)
        #macs
        module load MACS
        macs2 callpeak --gsize 1.87e9 --nomodel -t out/${name}/${name}_sorted.bam -n out/${name}/${name} --nolambda --slocal 10000 -q .01
        ;&
    *)
        ;;
esac

cd $olddir


