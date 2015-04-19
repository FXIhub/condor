#!/usr/bin/env python
import os

folder = "../src"
apply_header(folder)

def apply_header(p):
    if os.path.isdir(p):
        print "Process directory: %s" % p
        for pp in os.listdir(p):
            apply_header(pp)
    elif p[-len(".py"):] == ".py":
        ll_new = ""
        with open(p,"r") as f:
            if ll[0][:len("#!")] == "#!":
                ll.pop(0)
            in_header = True    
            for i,l in zip(range(ll),ll):
                if l[0] != "#":
                    in_header = False
                if not in_header:
                    ll_new.append(l)
        with open(p,"w") as f:
            f.writelines(ll_new)
            print p


