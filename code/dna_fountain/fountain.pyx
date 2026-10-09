"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from droplet import Droplet
from math import ceil
from utils import screen_repeat
from lfsr import lfsr, lfsr32p, lfsr32s
from robust_solition import PRNG
from reedsolo import RSCodec
import operator
import sys

import random

class DNAFountain:

    def __init__(self, 
                file_in, 
                file_size, 
                chunk_size,
                alpha, 
                stop = None,
                rs = 0, 
                c_dist = 0.1, 
                delta = 0.5, 
                np = False,
                max_homopolymer = 3,
                gc = 0.05
                ):

        self.file_in = file_in
        self.chunk_size = chunk_size
        self.num_chunks = int(ceil(file_size / float(chunk_size)))
        self.file_size = file_size
        self.alpha = alpha
        self.stop = stop
        self.final = self.calc_stop()

        self.lfsr = lfsr(lfsr32s(), lfsr32p()) 
        self.lfsr_l = len(    '{0:b}'.format( lfsr32p() )   ) - 1 
        self.seed = next(self.lfsr)

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np) 
        self.PRNG.set_seed(self.seed)

        self.rs = rs 
        self.rs_obj = RSCodec(self.rs)

        self.gc = gc
        self.max_homopolymer = max_homopolymer
        self.tries = 0 
        self.good = 0 
        self.oligo_l = self.calc_oligo_length()

    def calc_oligo_length(self):
        
        bits = self.chunk_size * 8 + self.lfsr_l + self.rs * 8
        return bits//4

    def calc_stop(self):

        if self.stop is not None:
            return self.stop
        
        stop = int(self.num_chunks*(1+self.alpha))+1
        return stop

    def droplet(self):
        
        data = None

        d, num_chunks = self.rand_chunk_nums() 

        for num in num_chunks: 
            if data is None: 
                data = self.chunk(num) 
            else: 
                data = list(map(operator.xor, data, self.chunk(num)))

        self.tries +=  1 

        return Droplet(data = data, 
                       seed = self.seed, 
                       rs = self.rs,
                       rs_obj = self.rs_obj,
                       num_chunks = num_chunks,
                       degree = d)

    def chunk(self, num):
        
        return self.file_in[num]

    def updateSeed(self):
        
        self.seed = next(self.lfsr) 
        self.PRNG.set_seed(self.seed) 

    def rand_chunk_nums(self):

        self.updateSeed() 
        blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
        return d, ix_samples 

    def screen(self, droplet):

        if screen_repeat(droplet, self.max_homopolymer, self.gc):
        
            self.good += 1
            return 1
        return 0

