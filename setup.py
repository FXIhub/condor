#!/usr/bin/env python
# -----------------------------------------------------------------------------------------------------
# CONDOR
# Simulator for diffractive single-particle imaging experiments with X-ray lasers
# http://xfel.icm.uu.se/condor/
# -----------------------------------------------------------------------------------------------------
# Copyright 2016 Max Hantke, Filipe R.N.C. Maia, Tomas Ekeberg
# Condor is distributed under the terms of the BSD 2-Clause License
# -----------------------------------------------------------------------------------------------------
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------------------------------
# General note:
# All variables are in SI units by default. Exceptions explicit by variable name.
# -----------------------------------------------------------------------------------------------------

# setup.py (simplified)

from setuptools import setup, Extension
import os
import numpy

# --- C Extension Build Logic ---

# Enable using threads (requires nfft installation with threads)
# This is now controlled ONLY by environment variables.
ENABLE_THREADS = bool(os.environ.get("CONDOR_ENABLE_THREADS"))

# Specify the include/library directories of the NFFT library via environment variables
NFFT_LIBRARY_DIR = os.environ.get("NFFT_LIBRARY_DIR")
NFFT_INCLUDE_DIR = os.environ.get("NFFT_INCLUDE_DIR")

print(numpy.get_include())
def get_extensions():
    """Constructs the list of C extensions to build."""

    # 1. Simple Icosahedron Extension
    ext_icosahedron = Extension(
        "condor.utils.icosahedron",
        sources=[os.path.join('condor', 'utils', 'icosahedronmodule.c')],
        include_dirs=[numpy.get_include()],
    )

    # 2. Conditional NFFT Extension
    # This logic is preserved from the original setup.py but simplified
    # to rely only on environment variables.
    
    library_dirs = [NFFT_LIBRARY_DIR] if NFFT_LIBRARY_DIR else []
    include_dirs = [numpy.get_include()]
    if NFFT_INCLUDE_DIR:
        include_dirs.append(NFFT_INCLUDE_DIR)

    if ENABLE_THREADS:
        libraries = ["nfft3_threads", "fftw3_threads", "fftw3"]
        define_macros = [("ENABLE_THREADS", None)]
    else:
        libraries = ["nfft3"]
        define_macros = []

    # The -rpath argument helps the dynamic linker find the shared libraries at runtime
    extra_link_args = []
    if library_dirs:
        rpath_arg = "-Wl," + ",".join(f"-rpath,{d}" for d in library_dirs)
        extra_link_args.append(rpath_arg)
        
    ext_nfft = Extension(
        "condor.utils.nfft",
        sources=[os.path.join('condor', 'utils', 'nfftmodule.c')],
        library_dirs=library_dirs,
        libraries=libraries,
        include_dirs=include_dirs,
        define_macros=define_macros,
        runtime_library_dirs=library_dirs, # For some linkers
        extra_link_args=extra_link_args,
    )

    return [ext_icosahedron, ext_nfft]

# --- Setup Call ---
# The 'setup()' call is now very minimal.
# All metadata is in pyproject.toml. Setuptools reads it automatically.
# We only need to provide the extension modules here.
setup(
    ext_modules=get_extensions(),
)