from bioutilities import Coordinate
import numpy as np
from bx.intervals.intersection import Intersecter, Interval
import os, subprocess
from scipy.io.matlab import savemat

wig_path=''
wig_mask='.phastCons44way.hg19.compiled'
random_factor=100
input_bed_file='H3K27ME3_ALL_cells_without_chrX_all_positive_coords.bed'
exon_bed_file='exons_hg19.bed'
intron_bed_file='introns_hg19.bed'
intergenic_bed_file='intergenic_regions_hg19.bed'

def intersections_lengths(n_input_coordinates,interval_tree,coord_to_row_index,target_coordinates):
    inters_lenghts=np.zeros(n_input_coordinates)
    for cl_index,c in enumerate(target_coordinates):

        if interval_tree.has_key(c.chr_id):
            coords_hits=interval_tree[c.chr_id].find(c.bpstart, c.bpend)

            for coord_hit in coords_hits:
                c_to_add=Coordinate.coordinates_from_interval(c.chr_id, coord_hit)
                row_index=coord_to_row_index[c_to_add]
                inters_lenghts[row_index]+=len(c &  c_to_add)
                inters_lenghts[row_index]=min(inters_lenghts[row_index],len(c_to_add)) #shared regions across genes..
                #print 'target:'+str(c),'hit:'+str(c_to_add),len(c),len(c_to_add),len(c &  c_to_add)

    return inters_lenghts

def extract_random_subcoordinate(c,length):
    bpstart=np.random.randint(c.bpstart,c.bpend-length+2)
    bpend=bpstart+length-1
    return Coordinate(c.chr_id,bpstart,bpend)

def extract_random_scores(c,target_coordinates,random_factor):

    random_sampling_score=np.zeros(random_factor)
    
    for idx_rnd in range(random_factor):

        not_founded=True
        while not_founded:
            c_random=target_coordinates[np.random.randint(len(target_coordinates))]
            not_founded= len(c_random) < len(c)

        region_nok=True
        while region_nok:
            try:
                region_nok=False
                c_en_random=extract_random_subcoordinate(c_random,len(c))
                #print c_en_random, len(c_en_random), len(c)
                random_sampling_score[idx_rnd]=read_from_wig(c_en_random,wig_path,wig_mask=wig_mask,only_average=True)[1]
            except:
                #print 'problema coordinate random'
                region_nok=True

                not_founded=True
                while not_founded:
                    c_random=target_coordinates[np.random.randint(len(target_coordinates))]
                    not_founded= len(c_random) < len(c)
                

    #print random_sampling_score
    return random_sampling_score


print 'carica i dati'

input_coordinates=Coordinate.bed_to_coordinates(input_bed_file)
exons_coordinates=Coordinate.bed_to_coordinates(exon_bed_file)
introns_coordinates=Coordinate.bed_to_coordinates(intron_bed_file)
intergenic_coordinates=Coordinate.bed_to_coordinates(intergenic_bed_file)

print 'alloca memoria per lunghezza intersezioni'
inters_length_exon_intron_intergenic=np.zeros((len(input_coordinates),3))

print 'costruisci interval tree degli enanchers'
interval_tree=dict()
coord_to_row_index=dict()
row_index=0
for c in input_coordinates:
    if c.chr_id not in interval_tree:
        interval_tree[c.chr_id]=Intersecter()

    interval_tree[c.chr_id].add_interval(Interval(c.bpstart,c.bpend))
    coord_to_row_index[c]=row_index
    row_index+=1

#per ogni categoria controlla intersezioni e update la riga con max intersezione

inters_length_exon_intron_intergenic[:,0]=intersections_lengths(len(input_coordinates),interval_tree,coord_to_row_index,exons_coordinates)
inters_length_exon_intron_intergenic[:,1]=intersections_lengths(len(input_coordinates),interval_tree,coord_to_row_index,introns_coordinates)
inters_length_exon_intron_intergenic[:,2]=intersections_lengths(len(input_coordinates),interval_tree,coord_to_row_index,intergenic_coordinates)

#per ogni elemento in enancher prendi il massimo per decidere la categoria
labels=inters_length_exon_intron_intergenic.argmax(1)


input_scores=np.zeros(len(input_coordinates))
random_scores=np.zeros((len(input_coordinates),random_factor))  


for idx,c in enumerate(input_coordinates):
    print idx, c
    try:
        input_scores[idx]=read_from_wig(c,wig_path,wig_mask=wig_mask,only_average=True)[1]

        if labels[idx] ==0:
            random_scores[idx,:]=extract_random_scores(c,exons_coordinates,random_factor)

        elif labels[idx]==1:
            random_scores[idx,:]=extract_random_scores(c,introns_coordinates,random_factor)

        else:
            random_scores[idx,:]=extract_random_scores(c,intergenic_coordinates,random_factor)
    except:
        print 'problema in:',idx,c,input_scores[idx],labels[idx]
    

np.save(prefix_to_add+'_scores',input_scores)
np.save(prefix_to_add+'random_scores',random_scores)
np.save(prefix_to_add+'labels',labels)
savemat(prefix_to_add+'risultati_conservation',{prefix_to_add+'labels':labels,prefix_to_add+'random_scores':random_scores,prefix_to_add+'_scores':input_scores})






