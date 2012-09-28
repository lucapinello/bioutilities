from bioutilities import Genome, Coordinate, Genome_mm
g=Genome('c:\data\genomes\hg18')
g1=Genome_mm('c:\data\genomes\hg18')

c=Coordinate('chr1',1,100)
print g.extract_sequence(c)

print g1.extract_sequence(c)

c=Coordinate('chr1',246249720,246249769,strand='-')

print g.extract_sequence(c)
print g1.extract_sequence(c)
