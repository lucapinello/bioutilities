import glob
import shlex
import re
import subprocess as sb
import shutil

names=set()
r=re.compile(r'^wgEncodeBroad(Histone|Tfbs)Gm12878(?P<name>.*)StdAln.*(1|2).bam')

for filename in glob.glob('*.bam'):
    print filename
    print r.match(filename).groupdict()['name']
    names.add(r.match(filename).groupdict()['name'])

print names

for name in names:
    print 'Processing:', name

    bam_file_1,bam_file_2=tuple(glob.glob('*'+name+'*.bam'))

    cmd='samtools view -H %s > current_header.sam' % bam_file_1
    process=sb.call(cmd,shell=True)
    print 'generated header'
    
    

    bam_filename=''.join([name,'_sample_merged_sorted'])
    cmd='samtools merge -h current_header.sam  -f /dev/stdout %s %s | samtools sort /dev/stdin %s -m 5000000000 ' %( bam_file_1,bam_file_2,bam_filename)
    print cmd
    process=sb.call(cmd,shell=True)
    print 'Processed:', name

print 'All Reps converted to bed'


