#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 18:07:12 2024

@author: joris
"""

import h5py
import numpy as np

n=100
z=300

I=np.random.rand(n,n,z)
File=h5py.File('essai.hdf5','w')

File.create_dataset('I',data=I,chunks=(50,50,50))
File.close()