import os
import glob


for filename in glob.glob('*.wig'):
    print filename
    os.system( ' '.join(['java -Xmx1500m  -jar igvtools.jar tile ',filename,filename.replace('.wig','tdf'),'hg19.genome']))
