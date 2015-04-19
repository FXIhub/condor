#!/usr/bin/env python
import os

HEADER = ["# -----------------------------------------------------------------------------------------------------\n",
          "# CONDOR\n",
          "# Simulator for diffractive single-particle imaging experiments with X-ray lasers\n",
          "# http://xfel.icm.uu.se/condor/\n",
          "# -----------------------------------------------------------------------------------------------------\n",
          "# Copyright 2014 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg\n",
          "# Condor is distributed under the terms of the GNU General Public License\n",
          "# -----------------------------------------------------------------------------------------------------\n",
          "# This program is free software; you can redistribute it and/or modify\n",
          "# it under the terms of the GNU General Public License as published by\n",
          "# the Free Software Foundation; either version 2 of the License, or\n",
          "# (at your option) any later version.\n",
          "# This program is distributed in the hope that it will be useful,\n",
          "# but without any warranty; without even the implied warranty of\n",
          "# merchantability or fitness for a pariticular purpose. See the\n",
          "# GNU General Public License for more details.\n",
          "# You should have received a copy of the GNU General Public License\n",
          "# along with this program; if not, write to the Free Software\n",
          "# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA\n",
          "# -----------------------------------------------------------------------------------------------------\n",
          "# General note:\n",
          "# All variables are in SI units by default. Exceptions explicit by variable name.\n",
          "# -----------------------------------------------------------------------------------------------------\n"]


str_starts_with = lambda s,s_start: s[:len(s_start)] == s_start
str_ends_with   = lambda s,s_end  : s[-len(s_end):]  == s_end

def apply_header(p):
    if os.path.isdir(p):
        print "Process directory: %s" % p
        for pp in os.listdir(p):
            apply_header(p+"/"+pp)
    elif str_ends_with(p,".py") and not str_starts_with(p.split("/")[-1],"."):
        ll_new = []
        with open(p,"r") as f:
            ll = f.readlines()
            if ll[0][:len("#!")] == "#!":
                ll_new.append(ll.pop(0))
            ll_new += HEADER
            in_header = True    
            for l in ll:
                if l[0] != "#":
                    in_header = False
                if not in_header:
                    ll_new.append(l)
        with open(p,"w") as f:
            f.writelines(ll_new)
            print p


folder = "../src"
apply_header(folder)
