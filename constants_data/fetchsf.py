import os, glob, numpy, pickle, sys, re

def generate_datafile(inpath,outpath): 
    path = './'

    sf_list = []


    PICKLEFILE = open('%s/elements.dat' % (outpath),'w')
    MASS_FILE = open('%s/masses.txt' % (inpath),'r')
    masslines = MASS_FILE.readlines()

    sys.stdout.write('[')
    dmasses = {}
    OUT = {}



    #for infile in glob.glob(os.path.join(path, 'sf/*.nff') ):
    for infile in glob.glob(os.path.join(path, '%s/*.nff' % (inpath)) ):
        #print "current file is: " + infile
        regexp = re.search("%s/([a-z]+).nff$" % (inpath),infile)
        el = regexp.group(1).capitalize()
        #el = infile[5:7].capitalize()
        OUT[el] = []
        FILE = open(infile)
        lines = FILE.readlines()
        lines = lines[1:]
        for line in lines:
            arg = line.split()
            if float(arg[1]) > 0:
                OUT[el].append([float(arg[0]),float(arg[1]),float(arg[2])])
        FILE.close()
        OUT[el] = numpy.array(OUT[el])

        symlen = len(el)
        for massline in masslines:
             if massline[4:4+symlen] == el:
                if massline[46] == '[':
                    if massline[49] == ']':
                        dmasses[el] = float(massline[47:49])
                    else:
                        dmasses[el] = float(massline[47:50])
                else:
                    dmasses[el] = float(massline[46:51])
                break

    pickle.dump([dmasses,OUT],PICKLEFILE)


    sys.stdout.write(']\n')



    MASS_FILE.close()
    PICKLEFILE.close()

if __name__ == "__main__":
    generate_datafiles("sf","..")
    
