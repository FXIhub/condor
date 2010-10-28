import spowpy
reload(spowpy)
error_counter = 0

totestlist = ["inp=spowpy.Input('spow.conf')","inp.set_printmode('stdout')","inp.set_printmode('loglist')","inp.set_plotmode('1d')","inp.set_plotmode('2d')","inp.set_sample_empty_densitymap(400E-09)","inp.sample.densitymap2d[1][1] = 1","if inp.sample.get_area(inp.sample.densitymap2d,inp.sample.densitymap_d) != inp.sample.densitymap_d**2: 1/0"]

try:
    print "(1) Running spow in default configuration simulating homogeneous sphere)."
    print "- create input object"
    inp = spowpy.Input()
    print "- run spow"
    out = spowpy.spow(inp)
    print "No errors occured."
except:
    print "!!! Error occured."
    error_counter += 1

try:
    print "(2) Running spow in default configuration simulating virus + goldball."
    print "- create input"
    inp = spowpy.Input()
    print "- create virus"
    inp.set_sample_virus_densitymap(100E-09,0.1,0.1,0.1)
    print "- put spheres"
    inp.sample.put_sphere(80E-09,200E-09,200E-09,cAu=1.0,massdensity=spowpy.DICT_massdensity['Au'])
    inp.sample.put_sphere(80E-09,400E-09,400E-09,materialtype='protein')
    print "- run spow"
    out = spowpy.spow(inp)
    print "No errors occured."
except:
    print "!!! Errors occured."
    error_counter += 1

print "(3) Testing all functions in 'totestlist'."
for command in totestlist:
    try:
        exec "%s" % command
        print "- %s leads to no errors." % command
    except:
        print "- !!! %s leads to error." % command
        error_counter += 1
print "During tests %i errors occured." % error_counter
        
