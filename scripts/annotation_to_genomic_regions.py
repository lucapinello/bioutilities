'''
To obtain the intergenic regions use hg19.genome (chr lengths) and then with bedtools
subtractBed -a hg19.bed -b coordinates_gene_hg19.bed > intergenic_regions_hg19.bed


'''

from bioutilities import Coordinate, Gene
gene_annotation_file='RefSeqhg19.txt'


gl=Gene.load_from_annotation(gene_annotation_file,load_exons_introns_info=True,header_lines=1)

exons=Gene.exons_from_annotations(gene_annotation_file)
Coordinate.coordinates_to_bed(exons,'hg_19_exons.bed',minimal_format=True)
del(exons)

introns=Gene.introns_from_annotations(gene_annotation_file)
Coordinate.coordinates_to_bed(introns,'hg_19_introns.bed',minimal_format=True)
del(introns)

Coordinate.coordinates_to_bed(Gene.genes_coordinates_from_annotations(gene_annotation_file),'genes_coordinates_hg_19.bed',minimal_format=True)
