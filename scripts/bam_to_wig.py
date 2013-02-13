import glob
import shlex, subprocess

for bam_filename in glob.glob('*.bam'):
    print bam_filename
    name=bam_filename.split('_')[0]+bam_filename.split('_')[1]
    #output_filename=bed_filename.replace('.bed','_merged.bed')
    cmd ='java  -Xmx2g  -jar  ~/home/programs/IGVTools/igvtools.jar  count  -w 100 -e 200 %s %s.wig  hg19' % (bam_filename,name)
    args = shlex.split(cmd)
    
    subprocess.Popen(cmd,shell=True)
    
