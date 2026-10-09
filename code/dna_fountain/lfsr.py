"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""

def lfsr(state, mask):
    
    result = state
    nbits = mask.bit_length()-1
    while True:
        result = (result << 1)
        xor = result >> nbits
        if xor != 0:
            result ^= mask

        yield result

def lfsr32p():

    return 0b100000000000000000000000011000101

def lfsr32s():

    return 0b001010101

def test():
    
    for pattern in lfsr(0b001, 0b100000000000000000000000011000101):
        print(pattern)
