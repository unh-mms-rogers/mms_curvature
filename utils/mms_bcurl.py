# This file adapted from mms_curl.py from the pyspedas library,
# sourced from https://github.com/spedas/pyspedas
#
# All modifications copyright 2020 Tim Rogers.  All rights reserved.
# Released under the MIT license.

import numpy as np
#from pytplot import get_data, store_data, options
#from pyspedas import tinterpol

#####
#Example of basic use within this package/library:
#
#from dataload_parallel import DataLoad
#
#numBirds = 4
#
#data_rate = 'srvy'
#trange = ['2017-06-17/20:00', '2017-06-17/21:00']
#dataset = DataLoad(trange=trange, data_rate=data_rate)
#fieldList = [None]*numBirds
#positionList = [None]*numBirds
#for bird in range(numBirds):
#    fieldList[bird] = (dataset['data']['fgm']['mms'+str(bird+1)+'_fgm_b_gsm_'+data_rate+'_l2'])
#    positionList[bird] = (dataset['data']['mec']['mms'+str(bird+1)+'_mec_r_gsm'])
#curl_dataset = mms_bcurl(fields=fieldList, positions=positionList)
#

#def mms_bcurl(fields=None, positions=None, suffix='', norm=False):
def mms_bcurl(postimes=None, posvalues=None, magtimes=None, magvalues=None, normalize=False):
    
    '''
    Calculates spacial gradient and curvature vector of the magnetic field.  
    Returns those and the master time (interpolated from FGM data) as numpy arrays

    Input parameters-
    NOTE:   order of each list is assumed to be consistent, i.e. that MMS1 data 
            is always in the [0] index of each list while MMS4 data is always 
            in the [3] index and so on.  While order is arbitrary, it must be 
            consistent for all of the following lists.

    postimes:   list of numpy arrays with the unix timestamp of each position vector

    posvalues:  list of numpy arrays with the position vector of a spacecraft

    magtimes:   list of numpy arrays with unix timestamp of each magnetic field vector

    magvalues:  list of numpy arrays with magnetic field vector measured by a
                spacecraft

    normalize:  if True normalizes magnetic field vectors before continusing
                calculation; required for calculating curvature

                if False leaves magnetic field as full vector with magnitude;
                required for calculating curl and divergence
    Notes:
        The input B-field data and position data are required to be in 
        orthogonal Cartesian coordinates such as GSE or GSM and must be the 
        coordinate system for all inputs.  Output will be in same coordinate 
        system as inputs.
 
        Based on the original mms_curl, written in IDL, by Jonathan Eastwood and 
        Malcolm Dunlop

        For more info on this method, see:
          Chanteur, G., Spatial Interpolation for Four Spacecraft: Theory, 
          Chapter 14 of Analysis methods for multi-spacecraft data, G. 
          Paschmann and P. W. Daly (Eds.) ISSI Scientific Report SR-001. 

    Returns:
        Dictionary of variables calculated
    '''
    
    # Number of spacecraft in inputs.  Assumed from number of position time arrays.
    numBirds = len(postimes)

    # Sanity check.  Throw an error if number of time arrays and data arrays don't all match.
    if len(posvalues) != numBirds: raise ValueError('Number of position value arrays does not match number of position time arrays!')
    if len(magtimes) != numBirds: raise ValueError('Number of magnetic field time arrays does not match number of position time arrays!')
    if len(magvalues) != numBirds: raise ValueError('Number of magnetic field value arrays does not match number of position time arrays!')


    bn = [None]*numBirds
    if normalize:
        # normalize magnetic fields
        if magvalues[0].shape[-1] == 4:   # tests for |B| given as 4th vector in data
            for bird in range(numBirds):
                bn[bird] = magvalues[bird][:,0:3]/magvalues[bird][:,3,np.newaxis]
        else: 
            for bird in range(numBirds):  # else calculates |B| from the 3-vector given
                bn[bird] = magvalues[bird]/np.linalg.norm(magvalues[bird], axis=1).reshape(magvalues[bird].shape[0], 1)

    else:
        # use magnetic fields without normalizing
        for bird in range(numBirds):
            bn[bird] = magvalues[bird][:,0:3]

    # find probe with latest beginning point for magnetic field data
    firsttimes = []
    lasttimes = []
    for bird in range(numBirds):
        firsttimes.append(magtimes[bird][0])
        lasttimes.append(magtimes[bird][-1])
    mastersc = np.argmax(firsttimes)
    # find earliest end time for all space craft in range
    tend = min(lasttimes)
    tend_i = 0  # initialize counting index for finding ending index for master S/C
    
    # initialize master time sequence by trimming time from last S/C to start (mastersc)
    while magtimes[mastersc][tend_i] < tend: tend_i += 1
    t_master = magtimes[mastersc][0:(tend_i+1)]
    
    # master mag field data arr, with interpolated values
    # Magnetic field data, interpolated to the previously determined master time sequence
    barr=np.ndarray((numBirds,t_master.shape[0],3))
    for bird in range(numBirds):
        barr[bird,:,0] = np.interp(t_master, magtimes[bird], bn[bird][:,0])
        barr[bird,:,1] = np.interp(t_master, magtimes[bird], bn[bird][:,1])
        barr[bird,:,2] = np.interp(t_master, magtimes[bird], bn[bird][:,2])
    


    # Calculate average |B| for export
    Bvals = [None]*numBirds
    if magvalues[0].shape[-1] == 4:    # tests for |B| given as 4th vector in data
        for bird in range(numBirds):
            Bvals[bird] = np.interp(t_master, magtimes[bird], magvalues[bird][:,3])
    else:
        for bird in range(numBirds):   # else calculates |B| from the 3-vector given
            Bvals[bird] = np.interp(t_master, magtimes[bird], np.linalg.norm(magvalues[bird], axis=1))
    bmag = np.average(Bvals, axis=0)


    
    # master position data array, with interpolated value
    # Spacecraft position data, interpolated to the previously determined master time sequence
    rarr = np.ndarray((numBirds,t_master.shape[0],3))
    for bird in range(numBirds):
        rarr[bird,:,0] = np.interp(t_master, postimes[bird], posvalues[bird][:,0])
        rarr[bird,:,1] = np.interp(t_master, postimes[bird], posvalues[bird][:,1])
        rarr[bird,:,2] = np.interp(t_master, postimes[bird], posvalues[bird][:,2])
    
    # Now all magnetic fields and positional data are of the same cadence and at the same times for each index
    # Indices are: [s/c(0=mms1, 1=mms2, 2=mms3, 3=mms4), time_step, vector(0=x, 1=y, 2=z)]
    # ie.  rarr[<spacecraft>, <timestep_index>, <cartesian_component_of_vector>]
    # eg.  Y-position of mms4 at first time step:  rarr[3,0,1]



    # Calculate the positional offset of each bird, relative to mms1, for each timestep
    ## posoffset[i] == offset of mms(i+2), relative to mms1
    posoffset = np.zeros([3, len(t_master), 3])
    for i in range(3):
        posoffset[i] = rarr[i+1]-rarr[0]
    
    
    k = np.zeros([4,len(t_master),3])
    
    ## Calculates the barycentre reciprical vector constants k_a, as inferred from equation (14.7) of "Analysis Methods for Multi-Spacecraft Data"
    # k[bird] = T_(matrix_multiply( 1/array_of_scalar_denominators , T_(array_of_vector_numerators)))
    # Breakout of above (for k[1]):
    #    numerator = np.cross(posoffset[1], posoffset[2])
    #    partial_denom = np.matmul(posoffset[0], np.transpose(numerator))
    #    k[1] = np.transpose(np.multiply((1/np.diag(partial_denom)), np.transpose(numerator)))
    #
    # Yes, this looks like a mess.
    # It also produces the desired results without resorting to a python for loop.
    # Using an intermediate numerator variable to help save memory...
    ## Preserving the below version because it's moderately easier to read than the more efficient form.
    #k_num = np.transpose(np.cross(posoffset[1], posoffset[2]))
    #k[1] = np.transpose(np.multiply((1/np.diag(np.matmul(posoffset[0], k_num))), k_num))
    #k_num = np.transpose(np.cross(posoffset[0], posoffset[2]))
    #k[2] = np.transpose(np.multiply((1/np.diag(np.matmul(posoffset[1], k_num))), k_num))
    #k_num = np.transpose(np.cross(posoffset[0], posoffset[1]))
    #k[3] = np.transpose(np.multiply((1/np.diag(np.matmul(posoffset[2], k_num))), k_num))
    
    ## Much more processor and memory efficient method, using numpy's Einstein Summation function.
    k_num = np.cross(posoffset[1], posoffset[2])
    k[1] = np.transpose(np.multiply((1/np.einsum('ij,ij->i',posoffset[0], k_num)), np.transpose(k_num)))
    k_num = np.cross(posoffset[0], posoffset[2])
    k[2] = np.transpose(np.multiply((1/np.einsum('ij,ij->i',posoffset[1], k_num)), np.transpose(k_num)))
    k_num = np.cross(posoffset[0], posoffset[1])
    k[3] = np.transpose(np.multiply((1/np.einsum('ij,ij->i',posoffset[2], k_num)), np.transpose(k_num)))
    del k_num
    ## Per equation (14.10) of "Analysis Methods for Multi-Spacecraft Data", sum of all k_a = 0.  Skip the heavy calculations for our relative origin.
    k[0] = k[0]-(k[3]+k[2]+k[1])

    # Per equation (14.15) of "Analysis Methods for Multi-Spacecraft Data", estimated gradient should be sum by bird of product of barcentre reciprical and transpose(B-vector)
    # np.reshape is required to permit numpy to orient the matrices for multiplication.
    # Since we need to add this empty dimension anyway, we use it to perform the required transpose operation manually.
    gradB = np.add.reduce(np.matmul(
                                    k.reshape(k.shape+(1,)),
                                    barr.reshape(barr.shape[:-1]+(1, barr.shape[-1]))
                                    ))
    
    # Per reference implementation, curlB is the sum across birds of each bird's barycentre reciprical vector cross B-vector, for each timestep
    curlB = np.add.reduce(np.cross(k,barr))

    # Per reference implementation, divB is the sum across birds of the dot product for each bird's barycentre reciprical vector and B-vector, for each timestep.
    #divB = np.einsum('ij,ij->i',barr[0], k[0])+np.einsum('ij,ij->i',barr[1], k[1])+np.einsum('ij,ij->i',barr[2], k[2])+np.einsum('ij,ij->i',barr[3], k[3])
    ## Replacing this with a more efficient implementation.
    #The reshape at the end is largely cosmetic, as the output is otherwise shape: (timestep, 1, 1)
    divB_diag = np.matmul(
                                    k.reshape(k.shape[:-1]+(1, k.shape[-1])),
                                    barr.reshape(barr.shape+(1,))
                                    )
    divB = np.add.reduce(divB_diag).reshape((barr.shape[1],))

    ## Sanity checks
    LevCiv3 = np.zeros((3,3,3))
    LevCiv3[0,1,2] = LevCiv3[1,2,0] = LevCiv3[2,0,1] = 1
    LevCiv3[0,2,1] = LevCiv3[2,1,0] = LevCiv3[1,0,2] = -1
    assert np.allclose(divB,
                       np.trace(gradB, axis1=-2, axis2=-1),
                       atol=0), 'Calculated divergence differs from trace of calculated gradient!'
    assert np.allclose(curlB,
                       np.einsum('ijk,...jk', LevCiv3, gradB),
                       atol=0), 'Calculated curl differs from Levi-Civita permutated gradient!'
    ## End of sanity checks

    jvec = 1e-12*curlB/m0
    jmag = np.sqrt(np.square(jvec[:,0])+np.square(jvec[:,1])+np.square(jvec[:,2]))


    # create the output variables
    out_vars['t_master' + suffix ] = t_master
    out_vars['barcentre' + suffix ] = barycentre
    out_vars['gradB' + suffix ] = gradB
    out_vars['curlB' + suffix ] = curlB
    out_vars['divB' + suffix  ] = divB
    #out_vars['divB_grad' + suffix  ] = divB_grad
    out_vars['jvec' + suffix] = jvec
    out_vars['jmag' + suffix] = jmag
    out_vars['baryB' + suffix] = baryB

    # Preserving old 'options' below, mostly as a reference for units.  Just in case.
    #
    #options('baryb' + suffix, 'ytitle', 'baryb [nT]')
    #options('divB' + suffix, 'ytitle', 'div(B) [nT/km]')
    #options('curlB' + suffix, 'ytitle', 'curl(B) [nT/km]')
    #options('curlB' + suffix, 'Color', ['b', 'g', 'r'])
    #options('curlB' + suffix, 'legend_names', ['delBx', 'delBy', 'delBz'])
    #options('jtotal' + suffix, 'ytitle', 'J [A/m^2]')
    #options('jtotal' + suffix, 'Color', ['b', 'g', 'r'])
    #options('jtotal' + suffix, 'legend_names', ['Jx', 'Jy', 'Jz'])
    #options('jperp' + suffix, 'ytitle', 'Jperp [A/m^2]')
    #options('jperp' + suffix, 'Color', ['b', 'g', 'r'])
    #options('jperp' + suffix, 'legend_names', ['Jperpx', 'Jperpy', 'Jperpz'])
    #options('jpar' + suffix, 'ytitle', 'Jparallel [A/m^2]')

    return out_vars
