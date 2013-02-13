import glob
import subprocess as sb

for bam_filename in glob.glob('*.bam'):
    print bam_filename
    cmd='bamToBed -i %s > %s' % (bam_filename,bam_filename.replace('.bam','.bed'))
    sb.call(cmd,shell=True)
    print 'Processed:',bam_filename

print 'All Reps converted to bed'



