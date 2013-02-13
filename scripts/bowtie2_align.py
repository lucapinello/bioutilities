import os,sys
import subprocess as sb
import glob
import shlex

bowtie2_index='/cluster7/lpinello/Genomes/hg19/hg19'

print 'CHECK THAT YOU ARE USING THE RIGHT GENOME!'

for in_filename in glob.glob('*.fastq.gz'):
    name=in_filename.replace('.fastq.gz','')
    aligned_bam_sorted_name=in_filename.replace('.fastq.gz','_sorted')

    print 'PROCESSING:', name
    cmd='bowtie2 -t --local -x %s -p 8 -U %s | samtools view -bS - | samtools sort -m 10000000000 - %s' %(bowtie2_index,in_filename,aligned_bam_sorted_name)
    print 'Executing:',cmd
    process = sb.call(cmd,shell=True)
    print 'FINISHED TO CONVERT',name 


    
