import os, glob

suffix_to_add='hg19.compiled'

for score_fname in glob.glob('*.wigFix'):
    print score_fname
    fields=score_fname.split('.')
    chr_id=fields[0]
    score_type=fields[1]
   
    wig_filename='.'.join([chr_id,score_type,suffix_to_add,'wig'])
    wib_filename='.'.join([chr_id,score_type,suffix_to_add,'wib'])

    print wig_filename
    print wib_filename


    os.system(' '.join(['wigEncode',score_fname,wig_filename,wib_filename]))
