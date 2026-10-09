"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from utils import int_to_four, four_to_dna
import random
import struct
import numpy as np

class Droplet:
    def __init__(self, data, seed, num_chunks = None, rs = 0, rs_obj = None, degree = None):

        self.data = data
        self.seed = seed
        self.num_chunks = set(num_chunks)
        self.rs = rs
        self.rs_obj = rs_obj
        self.degree = degree

        self.DNA = None

    def chunkNums(self):
        return self.num_chunks

    def toDNA(self, flag = None):

        if self.DNA is not None:
            return self.DNA
        
        self.DNA = int_to_four(self._package())
        return self.DNA

    def to_human_readable_DNA(self):
        
        return four_to_dna(self.toDNA())
        
    def _package(self):

        seed_ord =  [ c for c in struct.pack("!I", self.seed) ]

        message = seed_ord + list(self.data)
        
        if self.rs > 0:
            message = self.rs_obj.encode(message) 

        return message

