import os, glob, numpy, pickle, sys
 
path = './'

sf_list = []


PICKLEFILE = open('elements.dat','w')
MASS_FILE = open('../masses.txt','r')
masslines = MASS_FILE.readlines()

sys.stdout.write('[')
dmasses = []
OUT = []


for infile in glob.glob(os.path.join(path, '*.nff') ):
    #print "current file is: " + infile
    el = infile[2:]
    el = el[:-4].capitalize()
    descrSF = 'SF_' + el
    descrM = 'M_' + el
    vars()[descrSF] = []
    sys.stdout.write(descrSF + ',')
    FILE = open(infile)
    lines = FILE.readlines()
    lines = lines[1:]
    for line in lines:
        arg = line.split()
        if float(arg[1]) > 0:
            vars()[descrSF].append(numpy.array([float(arg[0]),float(arg[1])]))
    FILE.close()
    OUT.append(vars()[descrSF])
 
    symlen = len(el)
    for massline in masslines:
         if massline[4:4+symlen] == el:
            if massline[46] == '[':
                if massline[49] == ']':
                    dmasses.append([el,float(massline[47:49])])
                else:
                    dmasses.append([el,float(massline[47:50])])
            else:
                dmasses.append([el,float(massline[46:51])])
            break

pickle.dump([dmasses,OUT],PICKLEFILE)


sys.stdout.write(']\n')



MASS_FILE.close()
PICKLEFILE.close()
    
