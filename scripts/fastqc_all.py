import glob
import subprocess as sb
import os

for in_filename in glob.glob('*.fastq.gz'):
    print 'Processing:',in_filename
    output_folder='FASTQC'+in_filename.replace('.fastq.gz','')
    os.mkdir(output_folder)
    cmd='fastqc %s -o %s' %(in_filename,output_folder)
    sb.call(cmd,shell=True)
    print 'PROCESSED:',in_filename
