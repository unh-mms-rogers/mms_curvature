# Copyright 2020-2022 Anthony Rogers.  All rights reserved.
# Released under the Apache 2.0 license.

# This file is retained only for historical reference.

import numpy as np

def sigmaAB(A, uA, B, uB, op='*'):
    '''
    Calculates the uncertainty of two parameters with uncorrolated 
    errors propogated through either multiplication or division.
    
    For f = A * B   or   f = A / B, find sigma_f

    Parameters:

        A, B    - array-like (assumed to be 1D arrays) parameters used in the 
                  function 'f'.  Must be same size or broadcastable.

        uA, uB  - array-like (assumed constants or 1D arrays) of the uncertainty
                  associated with each measurement of the parameter A or B. 
                  Must be same size as or broadcastable to associated parameter
                  (uA to A; uB to B)

        op      - 'Operation'  Either '*' or '/' for if 
                  f = A * B
                  or   f = A / B

    Returns:
        'uf' array-like (assumed to be 1D array) of the propogated uncertainty 
        of the function 'f' in units of the function 'f' such that 'F = f +/- uf'

    See https://en.wikipedia.org/wiki/Propagation_of_uncertainty for more resources

    (c) ajr - 06.05.2020
    Apache 2.0 - See github.com/unh-mms-rogers/mms-curvature
    '''

    uAB = np.sqrt(np.power(np.divide(uA,A), 2) + np.power(np.divide(uB,B),2))

    if op == '*': f = np.multiply(A,B)
    elif op == '/': f = np.divide(A,B)
    else: 
        print("Invalid operator.  Must be '*' or '/'")
        return None

    return uAB, np.abs(np.multiply(uAB, f))


# end
