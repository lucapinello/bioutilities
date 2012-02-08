'''
Download chr lenghts from genome browser for example:
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e  "select chrom, size from hg19.chromInfo" > hg19.genome


'''


from bioutilities import Coordinate

chr_len_file='hg19.genome'
chr_lengths=dict()

with open(chr_len_file) as infile:

    infile.readline()
    for line in infile:
        fields=line.split()
        chr_lengths[fields[0]]=int(fields[1])


coordinates=[]

for chr_id in chr_lengths.keys():
    coordinates.append(Coordinate(chr_id,1,chr_lengths[chr_id]))

Coordinate.coordinates_to_bed(sorted(coordinates),'hg19.bed')
    
