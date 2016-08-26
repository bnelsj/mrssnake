import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.uint32

ctypedef np.uint32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def create_depth_array(unsigned int [:, :] count, unsigned int [:, :] depth, unsigned int rlen=36):

    for i in range(count.shape[0] - rlen):
        for j in range(count.shape[1]):
            if count[i, j] != 0:
                for k in range(rlen):
                    depth[i+k, j] += count[i, j]
