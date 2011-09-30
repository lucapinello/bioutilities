'''
Created on May 3, 2011

@author: lpinello
'''
import os, glob, string
import xlrd
from scipy.stats import rv_discrete
import mmap
import math

#from Bio import SeqIO
#from twobitreader import TwoBitFile

class constant_minus_one_dict (dict):
    def __missing__ (self, key):
        return -1

class constant_n_dict (dict):
    def __missing__ (self, key):
        return 'n'

nt2int=constant_minus_one_dict({'a':0,'c':1,'g':2,'t':3})      
int2nt=constant_n_dict({0:'a',1:'c',2:'g',3:'t'})
nt_complement=dict({'a':'t','c':'g','g':'c','t':'a'})

def read_sequence_from_fasta(fin, bpstart, bpend, line_length=50.0):
    bpstart=bpstart-1
    
    fin.seek(0)
    fin.readline()  #read the first line; the pointer is at the second line

    nbp = bpend - bpstart
    offset = int( bpstart + math.floor(bpstart/line_length)) #assuming each line contains 50 characters; add 1 offset per line
    
    if offset > 0:
        fin.seek(int(offset),1)

    seq = fin.read(nbp+int(math.floor(nbp/line_length))+1)
    seq = seq.replace('\n','')

    if len(seq) < nbp: 
        print 'Coordinate out of range:',bpstart,bpend
    
    return seq[0:nbp].lower()


class Coordinate:
    def __init__(self,chr_id,bpstart,bpend,name='ND',score=1.0,strand='ND'):
        self.chr_id=chr_id
        self.bpstart=bpstart
        self.bpend=bpend
        self.name=name
        self.score=score
        self.strand=strand


    def bpcenter(self):
        return (self.bpstart+self.bpend)/2

    def __eq__(self,other):
        return (self.chr_id==other.chr_id) & (self.bpstart==other.bpstart) & (self.bpend==other.bpend)
    
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __lt__(self, other):
        if chr_id2int[self.chr_id]<chr_id2int[other.chr_id]:
            return True 
        elif chr_id2int[self.chr_id]>chr_id2int[other.chr_id]:
            return False
        elif  self.bpstart<other.bpstart:
            return True
        elif self.bpstart>other.bpstart:
            return False
        else: 
            return False
    
    def __gt__(self,other):
        return not ( self.__lt__(other) or  self.__eq__(other)) 
    
    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)
    
    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)
            
    def __hash__(self):
        return hash((self.chr_id,self.bpstart,self.bpend))

    def __str__(self):
        return self.chr_id+'_'+str(self.bpstart)+'_'+str(self.bpend)+'_'+self.name+'_'+str(self.score)+'_'+self.strand
    
    def __repr__(self):
        return self.chr_id+'_'+str(self.bpstart)+'_'+str(self.bpend)+'_'+self.name+'_'+str(self.score)+'_'+self.strand
    
    def __len__(self):
        return self.bpend-self.bpstart+1
    
    
    @classmethod
    def read_coordinates_from_xls(cls,filename, chr_id_cl, bpstart_cl, bp_end_cl,name_cl=None,score_cl=None, strand_cl=None,header_lines=1):
        try:
            wb = xlrd.open_workbook(filename)
            sh = wb.sheet_by_index(0)
            chr_id_column = sh.col_values(chr_id_cl)[header_lines:]
            bpstart_column = sh.col_values(bpstart_cl)[header_lines:]
            bpend_column = sh.col_values(bp_end_cl)[header_lines:]
    
            coordinates = list()
            
            if name_cl:
                name_column = sh.col_values(name_cl)[header_lines:]
            else:
                name_column=['ND']*len(chr_id_column)
            
            if score_cl:
                score_column = sh.col_values(score_cl)[header_lines:]
            else:
                score_column=[1.0]*len(chr_id_column)
            
            if strand_cl:
                strand_column = sh.col_values(strand_cl)[header_lines:]
            else:
                strand_column=['ND']*len(chr_id_column)
                

            for chr_id, bpstart, bpend, name,score,strand in zip(chr_id_column, bpstart_column, bpend_column, name_column, score_column,strand_column):
                coordinates.append(Coordinate(str(chr_id), int(bpstart), int(bpend),name,float(score),strand))       
    
            return coordinates
        except:
            print "Error missing file or wrong columns number"
    
    @classmethod
    def bed_to_coordinates(cls,bed_file,header_lines=0, cl_chr_id=0, cl_bpstart=1, cl_bpend=2, cl_name=3, cl_score=4, cl_strand=5):
        with open(bed_file,'r') as infile:
            coordinates = list()
            
            for _ in range(header_lines):
                infile.readline()
            
            for line in infile:
                try:
                    coord=line.split()
                    chr_id=coord[cl_chr_id]
                    bpstart=coord[cl_bpstart]
                    bpend=coord[cl_bpend]
                    try:
                        name=coord[cl_name]
                    except:
                        name='ND'
                    try:
                        score=float(coord[cl_score])
                    except:
                        score=0                        
                    try:
                        strand=coord[cl_strand]
                    except:
                        strand='ND'
                    coordinates.append(Coordinate(str(chr_id), int(bpstart), int(bpend), name, score,strand))
                except:
                    print 'Skipping line:',line

            return coordinates
    
    @classmethod
    def bed_to_coordinates_dict(cls,bed_file,header_lines=0,cl_chr_id=0, cl_bpstart=1, cl_bpend=2, cl_name=3, cl_score=4, cl_strand=5):
        with open(bed_file,'r') as infile:
            coordinates = dict()
            
            for _ in range(header_lines):
                infile.readline()
            
            
            for idx,line in enumerate(infile):
                coord=line.split()
                chr_id=coord[cl_chr_id]
                bpstart=coord[cl_bpstart]
                bpend=coord[cl_bpend]
                try:
                    name=coord[cl_name]
                except:
                    name='ND'
                try:
                    score=coord[cl_score]
                except:
                    score='ND'                        
                try:
                    strand=coord[cl_strand]
                except:
                    strand='ND'
                c=Coordinate(str(chr_id), int(bpstart), int(bpend), name, score,strand)
                print c
                coordinates[c]=idx+1

            return coordinates
    
    @classmethod
    def coordinates_to_bed(cls,coordinates,bed_file):
        with open(bed_file,'w+') as outfile:
            for c in coordinates:
                outfile.write('%s\t%d\t%d\t%s\t%f\t%s\n' %(c.chr_id,c.bpstart,c.bpend,c.name,c.score,c.strand) )
    @classmethod
    def coordinates_to_nscore_format(cls,coordinates,bed_file):
        with open(bed_file,'w+') as outfile:
            for c in coordinates:
                outfile.write('%s\t%d\n' %(c.chr_id,c.bpcenter) )
    
    @classmethod
    def coordinates_from_interval(cls,chr_id,interval):
        return Coordinate(chr_id, interval.start,interval.end,)


    bpcenter=property(bpcenter)

class Gene:
    
    regions=[8000,2000,1000,1000]
    
    def __init__(self,access,name,coordinate,regions=None):
        self.access=access
        self.c=coordinate
        self.name=name
        
        if regions:
            self.regions=regions
    
    
    def __str__(self):
        return self.name+'_'+self.access+'_'+str(self.c)
    
    def tss(self):
        if self.c.strand=='-':
            return self.c.bpend
        else:
            return self.c.bpstart
    
    def end(self):
        if self.c.strand=='+':
            return self.c.bpend
        else:
            return self.c.bpstart
    
    def distal_c(self):
        if self.c.strand=='+':
            return Coordinate(self.c.chr_id,self.tss-self.regions[0],self.tss-self.regions[1]-1,strand='+')
        else:
            return Coordinate(self.c.chr_id,self.tss+self.regions[1]+1,self.tss+self.regions[0],strand='-') 
        
    def promoter_c(self):
        if self.c.strand=='+':
            return Coordinate(self.c.chr_id,self.tss-self.regions[1],self.tss+self.regions[2]-1,strand='+')
        else:
            return Coordinate(self.c.chr_id,self.tss-self.regions[2]+1,self.tss+self.regions[1],strand='-') 
    
    def intra_c(self):
        if self.c.strand=='+':
            return Coordinate(self.c.chr_id,self.tss+self.regions[2],self.end+self.regions[3],strand='+')
        else:
            return Coordinate(self.c.chr_id,self.end-self.regions[3],self.tss-self.regions[2],strand='-')
    
    def full_c(self):
        if self.c.strand=='+':
            return Coordinate(self.c.chr_id,self.tss-self.regions[0],self.end+self.regions[3],strand='+')
        else:
            return Coordinate(self.c.chr_id,self.end-self.regions[3],self.tss+self.regions[0],strand='-') 
        
    @classmethod    
    def load_from_annotation(cls,gene_annotation_file):
        genes_list=[]
        with open(gene_annotation_file,'r') as genes_file:
            
            for gene_line in genes_file:
                fields=gene_line.split('\t')
                chr=fields[2]
                bpstart=int(fields[4])
                bpend=int(fields[5])
                strand=fields[3]
                access=fields[1]
                name=fields[-4]
                c=Coordinate(chr,bpstart,bpend,strand=strand)
                genes_list.append(Gene(access,name,c))
                
        return genes_list
       
    tss=property(tss)
    end=property(end)
    distal_c=property(distal_c)
    promoter_c=property(promoter_c)
    intra_c=property(intra_c)
    full_c=property(full_c)
            
class Sequence:
    
    def __init__(self,seq=''):
        self.seq=string.lower(seq)

    def set_bg_model(self,ACGT_probabilities):
        self.bg_model= rv_discrete(name='bg', values=([0, 1, 2,3], ACGT_probabilities))
 
    def __str__(self):
        return self.seq
    
    def __len__(self):
        return len(self.seq)
   
    # seq_complement
    def _reverse_complement(self):
        return "".join([nt_complement[c] for c in self.seq[-1::-1]])

    reverse_complement=property(_reverse_complement)
    
    @classmethod
    def reverse_complement(self,seq):
        return "".join([nt_complement[c] for c in seq[-1::-1]])
    
    @classmethod
    def generate_random(cls,n,bgmodel=rv_discrete(name='bg', values=([0, 1, 2,3], [0.2955, 0.2045, 0.2045,0.2955]))):
        int_seq=cls.bg_model.rvs(size=n)
        return ''.join([int2nt[c] for c in int_seq])
    
class Seq_set(dict):
    def add_name(self,name):
        self.name=name
    def add_seq(self,seq,cord):
        self[cord]=seq

    def remove_seq_from_cord(self,cord):
        if cord in self.keys():
            del self[cord]
        else:
            print('sequence not found')

class Genome:
    def __init__(self,genome_directory,release='ND',verbose=False):
        self.chr=dict()
        self.genome_directory=genome_directory
        self.release=release
        self.chr_len=dict()
        
        for infile in glob.glob( os.path.join(genome_directory, '*.fa') ):
            try:

                chr_id=os.path.basename(infile).replace('.fa','')
                self.chr[chr_id] = open(infile,'r')
                
                self.chr_len[chr_id]=0
                
                self.chr[chr_id].readline()
                for line in self.chr[chr_id]:
                    self.chr_len[chr_id]+=len(line.strip())
                
                if verbose:    
                    print 'Read:'+ infile
            except:
                if verbose:
                    print 'Error, not loaded:',infile
        
        if verbose:
            print 'Genome initializated'
            
    def estimate_background(self):
        counting={'a':0,'c':0,'g':0,'t':0}

        for chr_id in self.chr.keys():
            if verbose:
                print 'Counting on:',chr_id
            

            self.chr[chr_id].seek(0)
            self.chr[chr_id].readline()
            
            for line in self.chr[chr_id]:
                for nt in counting.keys():
                    counting[nt]+=line.lower().count(nt)
        if verbose:
            print counting
        
        return counting

    
    def extract_sequence(self,coordinate, line_length=50.0):
        if not self.chr.has_key(coordinate.chr_id):
            if verbose:
                print "Warning: chromosome %s not present in the genome" % coordinate.chr_id
        else:

            bpstart=  coordinate.bpstart-1
            bpend=coordinate.bpend
            
            self.chr[coordinate.chr_id].seek(0)
            self.chr[coordinate.chr_id].readline()
            
           
            nbp = bpend - bpstart
            offset = int( bpstart + math.floor(bpstart/line_length))-1 

            if offset > 0:
                self.chr[coordinate.chr_id].seek(offset,1)
            
            seq = self.chr[coordinate.chr_id].read(nbp+int(math.floor(nbp/line_length))+1)
            seq = seq.replace('\n','')
            
            if len(seq) < nbp: 
                if verbose:
                    print 'Warning: coordinate out of range:',bpstart,bpend
            
            return seq[0:nbp].lower()            
            

class Genome_mm:
    
    def __init__(self,genome_directory,release='ND',verbose=False):
        self.chr=dict()
        self.genome_directory=genome_directory
        self.release=release        
        self.chr_len=dict()
        

        for infile in glob.glob( os.path.join(genome_directory, '*.fa') ):
            mm_filename=infile.replace('.fa','.mm')
            
            filename=infile.replace(genome_directory,'').replace('.fa','')
            chr_id=os.path.basename(infile).replace('.fa','')

            if not os.path.isfile(mm_filename):
                if verbose:
                    print 'Missing:'+chr_id+' generating memory mapped file (This is necessary only the first time) \n'
                with open(infile) as fi:
                    with open(os.path.join(genome_directory,chr_id+'.mm'),'w+') as fo:
                        #skip header
                        fi.readline()
                        for line in fi:
                            fo.write(line.strip()) 
                    if verbose:
                        print 'Memory mapped file generated for:',chr_id 

            with open(mm_filename,'r+') as f:
                self.chr[chr_id]= mmap.mmap(f.fileno(),0) 
                self.chr_len[chr_id]=len(self.chr[chr_id])
                if verbose:
                    print "Chromosome:%s Read"% chr_id
    
        if verbose:
            print 'Genome initializated'
        
        
    def extract_sequence(self,coordinate):
        return self.chr[coordinate.chr_id][coordinate.bpstart-1:coordinate.bpend].lower()
    
    
    def estimate_background(self):
        counting={'a':0,'c':0,'g':0,'t':0}

        for chr_id in self.chr.keys():
            if verbose:
                print 'Counting on:',chr_id
            counting[nt]+=self.chr[chr_id].lower().count(nt)
        
        if verbose:
            print counting
        
        return counting



''' OLD STUFF
def set_genome(genome='human'):
    if genome=='human':
        int2chr_id=dict(enumerate(['chr'+id for id in (map(str,range(1,23))+['X','Y'])],start=1));
    elif genome=='mouse':
        int2chr_id=dict(enumerate(['chr'+id for id in (map(str,range(1,20))+['X','Y','M'] ) ],start=1));
    else:
        raise Exception('not implemented')
    
    chr_id2int=dict((v,k) for k, v in int2chr_id.iteritems())
    
    return int2chr_id,chr_id2int

#HUMAN
#int2chr_id=dict(enumerate(['chr'+id for id in (map(str,range(1,23))+['X','Y'])],start=1));

#MOUSE
#int2chr_id=dict(enumerate(['chr'+id for id in (map(str,range(1,20))+['X','Y','M'] ) ],start=1));
#chr_id2int=dict((v,k) for k, v in int2chr_id.iteritems());


int2chr_id,chr_id2int=set_genome()
#print int2chr_id,chr_id2int

class Genome_2bit():
    
    #def __init__(self, 2bit_genome_file,release='ND'):
    def __init__(self,genome_file,release='ND'):
        self.data=TwoBitFile(genome_file)    
        self.release=release
    
    def extract_sequence(self,c):
        try:
            return self.data[c.chr_id][c.bpstart-1:c.bpend]
        except:
            print 'Bad coordinate in genome:',str(c)
        

class Genome:
    def __init__(self,genome_directory,number_of_chromosomes,release='ND'):
        self.chr=[None]*(number_of_chromosomes+1)
        self.genome_directory=genome_directory
        try:
            self.chr_len=map(int,open(os.path.join(genome_directory, 'chrlen.txt')).readlines())
            self.chr_len=[None]+self.chr_len
        except:
            self.chr_len=None
            

        self.release=release
        #print self.chr_len
        #print 'list_length:',len(self.chr_len)
        
        for infile in glob.glob( os.path.join(genome_directory, '*.fa') ):
            try:
                chr_idx=chr_id2int[ infile.replace(genome_directory,'').replace('.fa','')]
                self.chr[chr_idx] = open(infile,'r')
            except:
                print 'not loaded:',infile
            #print 'Readed:'+ infile
            #print 'length:',self.chr_len[chr_idx-1]
            
    def estimate_background(self):
        counting={'a':0,'c':0,'g':0,'t':0}

        for i in range(1,len(self.chr)):
            print self.chr[i]
            

            self.chr[i].seek(0)
            self.chr[i].readline()
            
            
            for line in self.chr[i]:
                for nt in counting.keys():
                    counting[nt]+=line.lower().count(nt)
        
        print counting
        return counting

                    
    #def __del__(self):
        #for file in chr:
        #    file.close()
       
    def read_seq_from_fasta(self,fin, bpstart, bpend):
        bpstart-=1
        bpend-=1
        fin.seek(0)
        head_length=len(fin.readline()) #assuming the first line is header
        nbp = int(bpend - bpstart + 1 + (bpend - bpstart + 1)/50) #devi considerare  pure gli a capo
        offset = bpstart + math.floor(bpstart/50) + head_length #assuming each line contains 50 characters; add 1 offset per line
        fin.seek(offset, 0)

        s = fin.read(nbp)
        
        if len(s) < nbp: 
            print 'Genome range out of scope for ',bpstart,bpend
   
        
    
    
    def read_sequence_from_fasta(self,fin, bpstart, bpend, line_length=50.0):
        bpstart=bpstart-1
        fin.seek(0)
        fin.readline()  #read the first line; the pointer is at the second line
        nbp = bpend - bpstart
        offset = bpstart + math.floor(bpstart/line_length) #assuming each line contains 50 characters; add 1 offset per line
        fin.seek(int(offset),1)
        seq = fin.read(nbp+int(math.floor(nbp/line_length))+1)
        seq = seq.replace('\n','')
        fin.close
        
        if len(seq) < nbp: 
            print 'Genome range out of scope for ',bpstart,bpend
        return seq.replace('\n','')[0:nbp].lower()
        
    
  
    def extract_sequence(self,coordinate):
        #if self.chr[chrstr2int(coordinate.chr_id)]==None:
        if self.chr[chr_id2int[coordinate.chr_id]]==None:
            print "Chromosome %s not readed"% coordinate.chr_id
        else:
            
            return self.read_seq_from_fasta(self.chr[chr_id2int[coordinate.chr_id]],coordinate.bpstart,coordinate.bpend)
            #return Sequence(self.chr[chrstr2int(coordinate.chr_id)][coordinate.bpstart-1:coordinate.bpend])
            #return Sequence(self.chr[chr_id2int[coordinate.chr_id]].seq[coordinate.bpstart-1:coordinate.bpend])
  
    def extract_sequence(self,coordinate):
        #if self.chr[chrstr2int(coordinate.chr_id)]==None:
        if self.chr[chr_id2int[coordinate.chr_id]]==None:
            print "Chromosome %s not readed"% coordinate.chr_id
        else:
            
            return self.read_sequence_from_fasta(self.chr[chr_id2int[coordinate.chr_id]],coordinate.bpstart,coordinate.bpend)
            #return Sequence(self.chr[chrstr2int(coordinate.chr_id)][coordinate.bpstart-1:coordinate.bpend])
            #return Sequence(self.chr[chr_id2int[coordinate.chr_id]].seq[coordinate.bpstart-1:coordinate.bpend])






class Genome_biopython:
    def __init__(self,genome_directory,number_of_chromosomes):
        self.chr=[None]*(number_of_chromosomes+1)
        self.genome_directory=genome_directory
        
        for infile in glob.glob( os.path.join(genome_directory, '*.fa') ):
            chr_idx=chr_id2int[ infile.replace(genome_directory,'').replace('.fa','')]
            print genome_directory
            print infile.replace(genome_directory,'').replace('.fa','')
            self.chr[chr_idx] = SeqIO.read(open(infile), "fasta")
            print 'Readed:'+ infile

class Genome_mmap:
        
    def extract_sequence(self,coordinate):
        #if self.chr[chrstr2int(coordinate.chr_id)]==None:
        if self.chr[chr_id2int[coordinate.chr_id]]==None:
            print "Chromosome %s not readed"% coordinate.chr_id
        else:
            #return Sequence(self.chr[chrstr2int(coordinate.chr_id)][coordinate.bpstart-1:coordinate.bpend])
            return Sequence(self.chr[chr_id2int[coordinate.chr_id]][coordinate.bpstart-1:coordinate.bpend])
            
     
    def fasta_to_mm(self,genome_directory):
        for infile in glob.glob( os.path.join(genome_directory, '*.fa') ):
            print "current file is: " + infile
            filename=infile.replace(genome_directory,'').replace('.fa','');
            with open(infile) as fi:
                with open(filename+'.mm','w+') as fo:
                    #skip header
                    fi.readline()
                    for line in fi:
                        fo.write(line.rstrip())  
    
    def read_chr(self,chr_idx):
        #file_to_read=genome_directory+'/'+int2chrstr(chr_idx)+'.mm'
        file_to_read=self.genome_directory+'/'+int2chr_id[chr_idx]+'.mm'
        with open(file_to_read,'r+') as f:
            self.chr[chr_idx]= mmap.mmap(f.fileno(),0) 
            print "Chromosome:%s Readed"% str(chr_idx)
    
    def read_all_chr(self):
        for chr_idx in range(1,len(self.chr)):
            #file_to_read=genome_directory+'/'+int2chrstr(chr_idx)+'.mm'
            file_to_read=self.genome_directory+'/'+int2chr_id[chr_idx]+'.mm'
            print file_to_read
            with open(file_to_read,'r+') as f:
                self.chr[chr_idx]= mmap.mmap(f.fileno(),0)
                print "Chromosome:%s Readed"% str(chr_idx)
        print "All Chromosome Readed"

'''
        
          
