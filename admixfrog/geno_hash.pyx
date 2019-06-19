"""trying to re-create admixtools hash functions
Not currently working
"""
from math import ceil
import itertools
from libc.stdio cimport *

def pyhash(l):                       
    cdef int h, hash_ = 0                         
    for idx in l:      
        h = 0                 
        for c in idx:                
            h = h * 23 + ord(c) #% 2 ** 32
        hash_ = hash_ * 17 #% 2 ** 32
        hash_ = hash_ ^ h #% 2 ** 32        
                                     
    return hash_                   



cdef expand_byte(byte b):
    cdef int b0, b1, b2, b3
    b0 = b // 64
    b1 = b % 64 // 16
    b2 = b % 16 // 4
    b3 = b % 4 
    return b0, b1, b2, b3

def read_geno(fname, max_snp=10000):
    cdef byte * row = NULL
    cdef FILE* f

    with open(fname, "rb") as f:
        head = f.read(48).strip(b"\x00").split()
        print(head)
        assert head[0] == b"GENO"
        n_ind, n_snp = int(head[1]), int(head[2])
        hash_ind, hash_snp = head[3:]

        rlen = row_length(n_ind)

        f.read(rlen - 48) # rest of header
        for i in range(min(n_snp, max_snp)):
            row = f.read(rlen)
            snp_iter = get_data(row, n_ind)
            yield itertools.islice(snp_iter, n_ind)

double cdef row_length(n_ind):
    return max(ceil(n_ind / 4 ), 48)

def get_data(row, n_ind):
    nb = ceil(n_ind / 4)
    for b in row[:nb]:
        for gt in expand_byte(b):
            yield gt
