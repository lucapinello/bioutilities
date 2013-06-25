from bioutilities import calculate_profile_matrix_bed_bam
import numpy as np
from scipy.io.matlab import savemat
import glob
import ntpath
import sys



table_filename=sys.argv[2]
'''
example of table format
H3K4me1	/gcdata/gcproj/Luca/Jennifer_GR/0/NS/BAM/JW10_RS4_11_No_Stim_H3K4me1_ChIP_010213.rmdup.bam	/gcdata/gcproj/Luca/Jennifer_GR/0/NS/BAM/JW7_RS4_11_No_Stim_MNase_Input_010213.rmdup.bam
H3K4me2	/gcdata/gcproj/Luca/Jennifer_GR/0/NS/BAM/JW13_RS4_11_No_Stim_H3K4me2_ChIP_010213.rmdup.bam	/gcdata/gcproj/Luca/Jennifer_GR/0/NS/BAM/JW7_RS4_11_No_Stim_MNase_Input_010213.rmdup.bam
.
.
.
'''

bed_filename=sys.argv[1]
try:
    WINDOW_SIZE=int(sys.argv[3])
except:
    WINDOW_SIZE=5000

try:
    if sys.argv[4]=='--use_input':
        USE_INPUT=True
    elif sys.argv[4]=='--no_input':
        USE_INPUT=False
    else:
        USE_INPUT=True
except:
    USE_INPUT=True
    


output_filename=bed_filename.replace('.bed','')+'IN'+table_filename.replace('.table','')
chip_profiles=dict()
chip_regions=dict()


chip_profiles_bg=dict()
chip_regions_bg=dict()


with open(table_filename) as table_file:

    for line in table_file:

        print 'Processing:',line.split('\t')[0]
        
        

        if USE_INPUT:

            chip_name,trg,inp=line.strip().split()
            print chip_name,trg,inp

            profile_matrix_trg,mapped_reads_trg=calculate_profile_matrix_bed_bam(bed_filename,trg,window_size=WINDOW_SIZE,use_strand=True)
            profile_matrix_inp,mapped_reads_inp=calculate_profile_matrix_bed_bam(bed_filename,inp,window_size=WINDOW_SIZE,use_strand=True)

            chip_profiles[chip_name]=    (profile_matrix_trg/float(mapped_reads_trg)).mean(0)  *1000000   
            chip_regions[chip_name]=    profile_matrix_trg/float(mapped_reads_trg)  *1000000  

            chip_profiles_bg[chip_name]= (profile_matrix_inp/float(mapped_reads_inp)).mean(0) *1000000
            chip_regions_bg[chip_name]= profile_matrix_inp/float(mapped_reads_inp) *1000000
 
        else:
            chip_name,trg=line.strip().split()
            print chip_name,trg


            profile_matrix_trg,mapped_reads_trg=calculate_profile_matrix_bed_bam(bed_filename,trg,window_size=WINDOW_SIZE, use_strand=True)
                        
            chip_profiles[chip_name]=(profile_matrix_trg/float(mapped_reads_trg)).mean(0)*1000000
            chip_regions[chip_name]=(profile_matrix_trg/float(mapped_reads_trg))*1000000

        print chip_profiles[chip_name]

    if USE_INPUT:
        savemat(output_filename,{'chip_profiles':chip_profiles,'chip_regions':chip_regions,'chip_profiles_bg':chip_profiles_bg,'chip_regions_bg':chip_regions_bg,'window_size':WINDOW_SIZE}) 
    else:
        savemat(output_filename,{'chip_profiles':chip_profiles,'chip_regions':chip_regions,'window_size':WINDOW_SIZE}) 
        
