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

def mms_bcurl(fields=None, positions=None, suffix=''):
    """
    This function applies the curlometer technique to MMS FGM data
    
    Parameters:
        fields : list of dict
            List of variables containing the B-field for each spacecraft 
            (in Cartesian coordinates, e.g. GSE or GSM)
            Each item in list is expected to be a data dictionary

        positions : list of dict
            List of variables containing the S/C position vectors for 
            each spacecraft (also in same coordinates as 'fields') 
            Each item in list is expected to be a data dictionary

        suffix: str
            The output variable names will be given this suffix.  By default, 
            no suffix is added.

    Notes:
        The input B-field data and position data are required to be in 
        orthogonal Cartesian coordinates such as GSE or GSM and must be the 
        coordinate system for all inputs.  Output will be in same coordinate 
        system as inputs.
 
        Based on the original mms_curl, written in IDL, by Jonathan Eastwood 

        For more info on this method, see:
          Chanteur, G., Spatial Interpolation for Four Spacecraft: Theory, 
          Chapter 14 of Analysis methods for multi-spacecraft data, G. 
          Paschmann and P. W. Daly (Eds.) ISSI Scientific Report SR-001. 

    Returns:
        Dictionary of variables calculated

    """

    if fields == None or positions == None:
        print('Error: B-field [fields] and spacecraft position [positions] keywords required.')
        return

    if len(fields) != 4 or len(positions) != 4:
        print('Error, fields and positions keywords should be specified as 4-element arrays containing the variable name for the field and position variables')
        return
    
    ## The following code to determine the timeseries to use was adapted from mms_curvature.py
    # find probe with latest beginning point for magnetic field data
    mastersc = np.argmax([fields[0]['x'][0], fields[1]['x'][0], fields[2]['x'][0], fields[3]['x'][0]])
    # find earliest end time for all space craft in range
    tend = min(fields[0]['x'][-1], fields[1]['x'][-1], fields[2]['x'][-1], fields[3]['x'][-1])
    tend_i = len(fields[mastersc]['x'])-1
    
    # initialize master time sequence by trimming time from last S/C to start (mastersc)
    while fields[mastersc]['x'][tend_i] > tend: tend_i -= 1
    timeseries = np.array(fields[mastersc]['x'][0:(tend_i+1)])
        
    out_vars = {}
    
    # *********************************************************
    # Magnetic Field
    # *********************************************************
    # interpolate the magnetic field data all onto the same timeline:
    # should be in GSE coordinates
    datab = np.ndarray((4,timeseries.shape[0],3))
    for bird in range(4):
        datab[bird,:,0] = np.interp(timeseries, fields[bird]['x'], fields[bird]['y'][:,0])
        datab[bird,:,1] = np.interp(timeseries, fields[bird]['x'], fields[bird]['y'][:,1])
        datab[bird,:,2] = np.interp(timeseries, fields[bird]['x'], fields[bird]['y'][:,2])

    # interpolate the definitive ephemeris onto the magnetic field timeseries
    # should be in GSE coordinates
    datapos = np.ndarray((4,timeseries.shape[0],3))
    for bird in range(4):
        datapos[bird,:,0] = np.interp(timeseries, positions[bird]['x'], positions[bird]['y'][:,0])
        datapos[bird,:,1] = np.interp(timeseries, positions[bird]['x'], positions[bird]['y'][:,1])
        datapos[bird,:,2] = np.interp(timeseries, positions[bird]['x'], positions[bird]['y'][:,2])

    m0 = 4.0*np.pi*1e-7 #  permeability of free space in SI units (H/m)

    # Calculate barycentre for each timestep
    barycentre = np.divide(np.add.reduce(datapos), datapos.shape[0])



    # Calculate the positional offset of each bird, relative to mms1, for each timestep
    ## posoffset[i] == offset of mms(i+2), relative to mms1
    posoffset = np.zeros([3, len(timeseries), 3])
    for i in range(3):
        posoffset[i] = datapos[i+1]-datapos[0]
    
    
    k = np.zeros([4,len(timeseries),3])
    
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
                                    datab.reshape(datab.shape[:-1]+(1, datab.shape[-1]))
                                    ))
    
    # Per reference implementation, curlB is the sum across birds of each bird's barycentre reciprical vector cross B-vector, for each timestep
    curlB = np.add.reduce(np.cross(k,datab))

    # Per reference implementation, divB is the sum across birds of the dot product for each bird's barycentre reciprical vector and B-vector, for each timestep.
    #divB = np.einsum('ij,ij->i',datab[0], k[0])+np.einsum('ij,ij->i',datab[1], k[1])+np.einsum('ij,ij->i',datab[2], k[2])+np.einsum('ij,ij->i',datab[3], k[3])
    ## Replacing this with a more efficient implementation.
    #The reshape at the end is largely cosmetic, as the output is otherwise shape: (timestep, 1, 1)
    divB_diag = np.matmul(
                                    k.reshape(k.shape[:-1]+(1, k.shape[-1])),
                                    datab.reshape(datab.shape+(1,))
                                    )
    divB = np.add.reduce(divB_diag).reshape((datab.shape[1],))

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
    out_vars['timeseries' + suffix ] = timeseries
    out_vars['barcentre' + suffix ] = barycentre
    out_vars['gradB' + suffix ] = gradB
    out_vars['curlB' + suffix ] = curlB
    out_vars['divB' + suffix  ] = divB
    out_vars['divB_grad' + suffix  ] = divB_grad
    out_vars['jvec' + suffix] = jvec
    out_vars['jmag' + suffix] = jmag

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
