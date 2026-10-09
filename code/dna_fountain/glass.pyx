from utils import *
from droplet import Droplet
from reedsolo import RSCodec
from robust_solition import PRNG
from other_screens import dexpandable_alphabet
import numpy as np
import operator
import sys
from collections import defaultdict
import pickle

class Glass:
    def __init__(self, num_chunks, out, header_size = 4, 
                 rs = 0, c_dist = 0.1, delta = 0.05, 
                flag_correct = True, gc = 0.2, max_homopolymer = 4, 
                max_hamming = 100, decode = True, chunk_size = 32, exDNA = False, np = False, truth = None):
        
        self.entries = []
        self.droplets = set()
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks
        self.header_size = header_size
        self.decode = decode
        self.chunk_size = chunk_size
        self.exDNA = exDNA
        self.np = np
        self.chunk_to_droplets = defaultdict(set)
        self.done_segments = set()
        self.truth = truth
        self.out = out

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np)

        self.max_homopolymer = max_homopolymer
        self.gc = gc
        prepare(self.max_homopolymer)
        self.max_hamming = max_hamming

        self.rs = rs
        self.RSCodec = None
        self.correct = flag_correct
        self.seen_seeds = set()

        if self.rs > 0:
            self.RSCodec = RSCodec(rs)

    def add_dna(self, dna_string):

        if self.exDNA:
            data = dexpandable_alphabet(dna_string, 
                                        len(dna_string), 
                                        n_symbols = 65,  
                                        n_bytes = 21, 
                                        alphabet_size = 6)
        else:
            data = dna_to_int_array(dna_string)

        if self.rs > 0:
            
            if self.correct: 
                try:
                    decode = self.RSCodec.decode(data)[0]
                    data_corrected = list(map(int, decode))

                except:
                    return -1, None 

                data_again = list(self.RSCodec.encode(data_corrected)) 

                if np.count_nonzero(data != list(data_again)) > self.max_hamming: 
                    
                    return -1, None

            else: 
                data_corrected  = data[0:len(data) - self.rs] 

        else:
            data_corrected = data

        seed_array = data_corrected[:self.header_size]
        seed = sum([   int(x)*256**i        for i, x in enumerate(seed_array[::-1])   ])
        payload = data_corrected[self.header_size:]

        if seed in self.seen_seeds:
            return -1, None
        self.add_seed(seed)

        if self.decode:

            self.PRNG.set_seed(seed)
            blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
            d = Droplet(payload, seed, ix_samples)

            if not screen_repeat(d, self.max_homopolymer, self.gc):
                return -1, None

            self.addDroplet(d)

        return seed, data

    def addDroplet(self, droplet):

        self.droplets.add(droplet)
        for chunk_num in droplet.num_chunks:
            self.chunk_to_droplets[chunk_num].add(droplet) 
        
        self.updateEntry(droplet) 

    def updateEntry(self, droplet):

        for chunk_num in (droplet.num_chunks & self.done_segments):

            droplet.data = list(map(operator.xor, droplet.data, self.chunks[chunk_num]))
            
            droplet.num_chunks.remove(chunk_num)
            
            self.chunk_to_droplets[chunk_num].discard(droplet)

        if len(droplet.num_chunks) == 1: 
            lone_chunk = droplet.num_chunks.pop() 
            self.chunks[lone_chunk] = droplet.data 
            self.done_segments.add(lone_chunk) 
            if self.truth:
                self.check_truth(droplet, lone_chunk)
            self.droplets.discard(droplet) 
            self.chunk_to_droplets[lone_chunk].discard(droplet) 

            for other_droplet in self.chunk_to_droplets[lone_chunk].copy():
                self.updateEntry(other_droplet)

    def getString(self):
        
        res = ''
        for x in self.chunks:
            res += ''.join(map(chr, x))
        return res

    def alive(self):
        return True

    def check_truth(self, droplet, chunk_num):
        try:
            truth_data = self.truth[chunk_num]
        except:
            print("Error. chunk:", chunk_num, " does not exist.")
            quit(1)

        if not droplet.data == truth_data:
            
            print("Decoding error in ", chunk_num, ".\nInput is:", truth_data,"\nOutput is:", droplet.data,"\nDNA:", droplet.to_human_readable_DNA(flag_exDNA = False))
            quit(1)
        else:
            
            return 1

    def add_seed(self, seed):
        self.seen_seeds.add(seed)

    def len_seen_seed(self):
        return len(self.seen_seeds)

    def isDone(self):
        if self.num_chunks - len(self.done_segments) > 0:
            return None 
        return True

    def chunksDone(self):
        return len(self.done_segments)

    def save(self):

        name =  self.out + '.glass.tmp'
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        return name

