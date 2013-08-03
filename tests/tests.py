import sys,os

print "Clean up"
os.system("cp ../setup.py .")
os.system("rm -r src")
os.system("cp -r ../src .")
os.system("rm -r const")
os.system("cp -r ../const .")
os.system("rm -r test_installation")
os.system("mkdir test_installation")

print "Install propagator"
here = os.path.dirname(os.path.realpath(__file__))
os.system("%s setup.py install --prefix=%s/test_installation" % (sys.executable,here))

print (here+"/test_installation/lib/python2.7/site-packages")
sys.path.insert(0,here+"/test_installation/lib/python2.7/site-packages")

import propagator
print propagator

I = propagator.Input()
O = propagator.propagator(I)
