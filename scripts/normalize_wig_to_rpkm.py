import glob

for input_file in glob.glob('*.wig'):

    with open(input_file) as infile:
        K=0.0
        for line in infile:
            try:
                K+=float(line.split()[1])
            except:
                print line
                pass

        print input_file,' normalization Factor:', K

        output_file=input_file.replace('.wig','_rpm_normalized.wig')
        infile.seek(0)
        with open(output_file,'w+') as outfile:
            for line in infile:
                    if not '#' in line:
                        fields=line.split()
                        if len(fields)==2:
                            normalized_value=(float(fields[1])/K)*1000000
                            outfile.write('%s\t%f\n' %( fields[0],normalized_value))
                        else:
                            print line
                            outfile.write(line)
                        
