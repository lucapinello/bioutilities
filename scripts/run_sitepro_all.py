import glob
import os
import subprocess


for wig_filename in  glob.glob('*.wig'):
    print wig_filename
    for bed_filename in glob.glob('*.bed'):
        print 'PROCESSING:',wig_filename, os.path.basename(bed_filename)
        
        name=wig_filename.replace('.wig','')+'_in_'+os.path.basename(bed_filename).replace('.bed','')
        print name
        cmd='sitepro -w %s -b %s --dump  --span 2500  --pf-res=100 --name %s' %(wig_filename, bed_filename, name )
        subprocess.check_call(cmd, shell=True)
        os.rename(wig_filename.replace('.wig','')+'_dump.txt',name+'_dump.txt')
        

   
