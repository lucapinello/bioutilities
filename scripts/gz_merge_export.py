import glob
import subprocess 
import shlex

experiments=set()
# scan to find the different experiment to merge
for infile in glob.glob('*.gz'):
    experiments.add(infile.split('-')[0])

print experiments

for exp in experiments:
    print 'Current Experiment:', exp
    command='zcat '
    for infile in glob.glob(exp+'*.gz'):
        print 'File to concat:',infile.split('_')[-2]
        command+=' %s' % infile
    command+=' > %s_export.txt && gzip %s_export.txt' %(exp,exp)
    print 'running:',command
    process = subprocess.Popen(shlex.split(command),shell=True)
    retcode = process.wait()
    print retcode
    
