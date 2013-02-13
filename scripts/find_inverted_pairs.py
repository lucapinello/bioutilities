import sys
infilename=sys.argv[1]
wrong_p_p=0
wrong_n_n=0
correct=0
first=False
n=0.0

with open(infilename) as infile:
    for line in infile:
        if '@' in line:
            pass
        else:
            fields=line.split('\t')
            r_id=fields[0]
            strand=(int(fields[1]) & 16 )>0
            
            if first:
                if r1_id==r_id:
		    n+=1
                    if r1_strand == strand:
                        if not r1_strand:
                            wrong_p_p+=1
                            print r_id+' + + on:',' '.join(fields[2:4])
                        else:
                            wrong_n_n+=1
                            print r_id+' - - on:',' '.join(fields[2:4])

                        
                        first=False
		    else:
			correct+=1
                else:
                    first=True
                    r1_id=r_id
                    r1_strand=strand

                
            else:
            
                r1_id=r_id
                r1_strand=strand
                first=True



print 'correct pairs %d %f' %(correct,correct/n)
print 'wrong pairs + + %d %f' %(wrong_p_p,wrong_p_p/n)
print 'wrong pairs - - %d %f' %(wrong_n_n,wrong_n_n/n)
print 'Hello Director!'
