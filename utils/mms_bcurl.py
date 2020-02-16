# This file adapted from mms_curl.py from the pyspedas library,
# sourced from https://github.com/spedas/pyspedas
#
# All modifications copyright 2020 Tim Rogers.  All rights reserved.
# Released under the MIT license.

import math
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
            (in GSE coordinates)
            Each item in list is expected to be a data dictionary

        positions : list of dict
            List of variables containing the S/C position vectors for 
            each spacecraft (also GSE coordinates) 
            Each item in list is expected to be a data dictionary

        suffix: str
            The output variable names will be given this suffix.  By default, 
            no suffix is added.

    Notes:
        The input B-field data and position data are required to be in 
        GSE coordinates
 
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
    timeseries = fields[mastersc]['x'][0:(tend_i+1)]
        
    out_vars = {}
    
    # *********************************************************
    # Magnetic Field
    # *********************************************************
    # interpolate the magnetic field data all onto the same timeline (MMS1):
    # should be in GSE coordinates
    ##tinterpol(fields[1], fields[0], newname=fields[1] + '_i')
    ##tinterpol(fields[2], fields[0], newname=fields[2] + '_i')
    ##tinterpol(fields[3], fields[0], newname=fields[3] + '_i')
    datab = np.ndarray((4,timeseries.shape[0],3))
    for bird in range(4):
        datab[bird,:,0] = np.interp(timeseries, fields[bird]['x'], fields[bird]['y'][:,0])
        datab[bird,:,1] = np.interp(timeseries, fields[bird]['x'], fields[bird]['y'][:,1])
        datab[bird,:,2] = np.interp(timeseries, fields[bird]['x'], fields[bird]['y'][:,2])
    #out_vars[fields[0] + '_i'] = np.interp(timeseries, fields[0]['x'], fields[0]['y'])
    #out_vars[fields[1] + '_i'] = np.interp(timeseries, fields[1]['x'], fields[1]['y'])
    #out_vars[fields[2] + '_i'] = np.interp(timeseries, fields[2]['x'], fields[2]['y'])
    #out_vars[fields[3] + '_i'] = np.interp(timeseries, fields[3]['x'], fields[3]['y'])

    # interpolate the definitive ephemeris onto the magnetic field timeseries
    # should be in GSE coordinates
    ##tinterpol(positions[0], fields[0], newname=positions[0] + '_i')
    ##tinterpol(positions[1], fields[0], newname=positions[1] + '_i')
    ##tinterpol(positions[2], fields[0], newname=positions[2] + '_i')
    ##tinterpol(positions[3], fields[0], newname=positions[3] + '_i')
    datapos = np.ndarray((4,timeseries.shape[0],3))
    for bird in range(4):
        datapos[bird,:,0] = np.interp(timeseries, positions[bird]['x'], positions[bird]['y'][:,0])
        datapos[bird,:,1] = np.interp(timeseries, positions[bird]['x'], positions[bird]['y'][:,1])
        datapos[bird,:,2] = np.interp(timeseries, positions[bird]['x'], positions[bird]['y'][:,2])
    #out_vars[positions[0] + '_i'] = np.interp(timeseries, positions[0]['x'], positions[0]['y'])
    #out_vars[positions[1] + '_i'] = np.interp(timeseries, positions[1]['x'], positions[1]['y'])
    #out_vars[positions[2] + '_i'] = np.interp(timeseries, positions[2]['x'], positions[2]['y'])
    #out_vars[positions[3] + '_i'] = np.interp(timeseries, positions[3]['x'], positions[3]['y'])

    m0 = 4.0*math.pi*1e-7

    ##timesb1, datab1 = get_data(fields[0])
    ##timesb2, datab2 = get_data(fields[1] + '_i')
    ##timesb3, datab3 = get_data(fields[2] + '_i')
    ##timesb4, datab4 = get_data(fields[3] + '_i')
    #datab1 = out_vars[fields[0] + '_i']
    #datab2 = out_vars[fields[1] + '_i']
    #datab3 = out_vars[fields[2] + '_i']
    #datab4 = out_vars[fields[3] + '_i']

    # extract the vector
    #b1 = datab[0][:, 0:3]
    #b2 = datab[1][:, 0:3]
    #b3 = datab[2][:, 0:3]
    #b4 = datab[3][:, 0:3]

    ##timesp1, p1 = get_data(positions[0] + '_i')
    ##timesp2, p2 = get_data(positions[1] + '_i')
    ##timesp3, p3 = get_data(positions[2] + '_i')
    ##timesp4, p4 = get_data(positions[3] + '_i')
    #p1 = datapos[0]
    #p2 = datapos[1]
    #p3 = datapos[2]
    #p4 = datapos[3]
    #p1 = out_vars[positions[0] + '_i']
    #p2 = out_vars[positions[1] + '_i']
    #p3 = out_vars[positions[2] + '_i']
    #p4 = out_vars[positions[3] + '_i']

    #divB = np.zeros([len(timeseries)])
    #curlB = np.zeros([len(timeseries), 3])
    #baryb = np.zeros([len(timeseries), 3])
    #baryb2 = np.zeros([len(timeseries), 3])
    #baryb3 = np.zeros([len(timeseries), 3])
    #baryb4 = np.zeros([len(timeseries), 3])
    #sampleb = np.zeros([len(timeseries), 3])

    #jtotal = np.zeros([len(timeseries), 4])
    #jtotal = np.zeros([len(timeseries), 3])
    #btotal = np.zeros([len(timeseries), 1])
    #jparallel = np.zeros([len(timeseries), 1])
    #jperpvec = np.zeros([len(timeseries), 4])
    #jperp = np.zeros([len(timeseries), 1])
    #alphaparallel = np.zeros([len(timeseries), 1])
    #alpha = np.zeros([len(timeseries), 1])

    ## posoffset[i] == offset of mms(i+2), relative to mms1
    posoffset = np.zeros([3, len(timeseries), 3])
    for i in range(3):
        posoffset[i] = datapos[i+1]-datapos[0]
    
    
    k = np.zeros([4,len(timeseries),3])
    
    ## Don't do this.  This gave bad results and occasionally crashed python...
    #k[1] = (np.transpose(np.cross(posoffset[1],posoffset[2]))*(1/np.matmul(posoffset[0], np.transpose(np.cross(posoffset[1],posoffset[2]))))).transpose()
    #k[2] = (np.transpose(np.cross(posoffset[0],posoffset[2]))*(1/np.matmul(posoffset[1], np.transpose(np.cross(posoffset[0],posoffset[2]))))).transpose()
    #k[3] = (np.transpose(np.cross(posoffset[0],posoffset[1]))*(1/np.matmul(posoffset[2], np.transpose(np.cross(posoffset[0],posoffset[1]))))).transpose()
    
    ## Calculates the barcentre reciprical vector constants k_a, as inferred from equation (14.7) of "Analysis Methods for Multi-Spacecraft Data"
    # k[bird] = T_(matrix_multiply( 1/array_of_scalar_denominators , T_(array_of_vector_numerators)))
    # Breakout of above (for k[1]):
    #    numerator = np.cross(posoffset[1], posoffset[2])
    #    partial_denom = np.matmul(posoffset[0], np.transpose(numerator))
    #    k[1] = np.transpose(np.multiply((1/np.diag(partial_denom)), np.transpose(numerator)))
    #
    #  Note to self:  numerator is only ever used by this method while transposed
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

    curlB = np.add.reduce(np.cross(k,datab))
    ## Replacing this with a more efficient implementation.
    #divB = np.diag(np.matmul(datab[0], np.transpose(k[0]))+np.matmul(datab[1], np.transpose(k[1]))+np.matmul(datab[2], np.transpose(k[2]))+np.matmul(datab[3], np.transpose(k[3])))
    divB = np.einsum('ij,ij->i',datab[0], k[0])+np.einsum('ij,ij->i',datab[1], k[1])+np.einsum('ij,ij->i',datab[2], k[2])+np.einsum('ij,ij->i',datab[3], k[3])

    jvec = 1e-12*curlB/m0
    jmag = np.sqrt(np.square(jvec[:,0])+np.square(jvec[:,1])+np.square(jvec[:,2]))


    # leave as a loop for now because you have to construct and manipulate a matrix for each time step.
    ### Tim:  Nevermind.  Don't need the loop anymore for the reduced calculations we're doing.
    ###       May need a loop back if we reinstate some of the other outputs.  Will need to look at vectorization of those if/when needed.
    #for i, time in enumerate(timeseries):
    #    # position of birds, relative to mms1
    #    p12 = p2[i, 0:3]-p1[i, 0:3]
    #    p13 = p3[i, 0:3]-p1[i, 0:3]
    #    p14 = p4[i, 0:3]-p1[i, 0:3]
    #
    #    
    #    # Calculates the barcentre reciprical vector constants k_a, as inferred from equation (14.7) of "Analysis Methods for Multi-Spacecraft Data"
    #    k2 = np.cross(p13, p14)*(1/(np.matmul(p12, np.transpose(np.cross(p13, p14)))))
    #    k3 = np.cross(p12, p14)*(1/(np.matmul(p13, np.transpose(np.cross(p12, p14)))))
    #    k4 = np.cross(p12, p13)*(1/(np.matmul(p14, np.transpose(np.cross(p12, p13)))))
    #
    #    # Per equation (14.10) of "Analysis Methods for Multi-Spacecraft Data", sum of all k_a = 0.  Skip the heavy calculations for our relative origin.
    #    k1 = 0-(k4+k3+k2)
    #    ###############
    #
    #    curlB[i] = np.cross(k1, b1[i, :])+np.cross(k2, b2[i, :])+np.cross(k3, b3[i, :])+np.cross(k4, b4[i, :])
    #    divB[i] = np.matmul(b1[i, :], k1) + np.matmul(b2[i, :], k2) + np.matmul(b3[i, :], k3) + np.matmul(b4[i, :], k4)
    #
    #    #gradbx = b1[i, 0]*k1 + b2[i, 0]*k2 + b3[i, 0]*k3 + b4[i, 0]*k4
    #    #gradby = b1[i, 1]*k1 + b2[i, 1]*k2 + b3[i, 1]*k3 + b4[i, 1]*k4
    #    #gradbz = b1[i, 2]*k1 + b2[i, 2]*k2 + b3[i, 2]*k3 + b4[i, 2]*k4
    #
    #    #barycentre = (p1[i, 0:3] + p2[i, 0:3] + p3[i, 0:3] + p4[i, 0:3])/4.0
    #
    #    ## and here is the field at the barycentre (calculate 4 ways)
    #    #baryb[i, 0] = b1[i, 0] + np.sum(gradbx*(barycentre-p1[i, 0:3]))
    #    #baryb[i, 1] = b1[i, 1] + np.sum(gradby*(barycentre-p1[i, 0:3]))
    #    #baryb[i, 2] = b1[i, 2] + np.sum(gradbz*(barycentre-p1[i, 0:3]))
    #
    #    #baryb2[i, 0] = b2[i, 0] + np.sum(gradbx*(barycentre-p2[i, 0:3]))
    #    #baryb2[i, 1] = b2[i, 1] + np.sum(gradby*(barycentre-p2[i, 0:3]))
    #    #baryb2[i, 2] = b2[i, 2] + np.sum(gradbz*(barycentre-p2[i, 0:3]))
    #
    #    #baryb3[i, 0] = b3[i, 0] + np.sum(gradbx*(barycentre-p3[i, 0:3]))
    #    #baryb3[i, 1] = b3[i, 1] + np.sum(gradby*(barycentre-p3[i, 0:3]))
    #    #baryb3[i, 2] = b3[i, 2] + np.sum(gradbz*(barycentre-p3[i, 0:3]))
    #
    #    #baryb4[i, 0] = b4[i, 0] + np.sum(gradbx*(barycentre-p4[i, 0:3]))
    #    #baryb4[i, 1] = b4[i, 1] + np.sum(gradby*(barycentre-p4[i, 0:3]))
    #    #baryb4[i, 2] = b4[i, 2] + np.sum(gradbz*(barycentre-p4[i, 0:3]))
    #
    #    # (these above all agree so this is the magnetic field at the barycentre)
    #    #divb[i, 0] = time
    #    #divb[i, 1] = divergence
    #    #divb[i, 2] = curlmag[0]
    #    #divb[i, 3] = curlmag[1]
    #    #divb[i, 4] = curlmag[2]
    #
    #    # the cross product of the calculated curl and the sample field times 1e-21 (SI), divided by m0
    #
    #    # curl is in nT/km, nT/km*1e-12 = T/m
    #    # field is in nT, nT*1e-9 = T
    #    # j is curl B / m0 (curl B = m0*j)
    #    # use the magnetic field at the barycentre
    #
    #    # compute the current components and total specifically
    #    #jtotal[i, 0:3] = 1e-12*divb[i, 2:5]/m0
    #    jtotal[i] = 1e-12*curlB[i]/m0
    #    #jtotal[i, 3] = np.sqrt(jtotal[i, 0]**2+jtotal[i, 1]**2+jtotal[i, 2]**2)
    #
    #    ## compute the parallel and perpendicular components of the current
    #    #btotal = np.sqrt(np.dot(baryb[i, 0:3], baryb[i, 0:3]))
    #
    #    ## parallel is J.B/|B|
    #    #jparallel[i] = np.dot(jtotal[i, 0:3], baryb[i, 0:3])/btotal
    #    #jparallel[i] = jparallel[i][0]
    #
    #    ## perp is J - J// B/|B| (components and total perpendicular current)
    #    #jperpvec[i, 0:3] = jtotal[i, 0:3] - (jparallel[i]*baryb[i, 0:3])/btotal
    #    #jperpvec[i, 3] = np.sqrt(jperpvec[i, 0]**2 + jperpvec[i, 1]**2 + jperpvec[i, 2]**2)
    #
    #    ## alpha parameter
    #    #alphaparallel[i] = np.abs(jparallel[i])/(1e-9*btotal)
    #    #alpha[i] = np.abs(jtotal[i, 3])/(1e-9*btotal)


    # create the output variables
    ##store_data('baryb' + suffix, data={'x': timesb1, 'y': baryb})
    ##store_data('curlB' + suffix, data={'x': timesb1, 'y': divb[:, 2:5]})
    ##store_data('divB' + suffix, data={'x': timesb1, 'y': divb[:, 1]})
    ##store_data('jtotal' + suffix, data={'x': timesb1, 'y': jtotal[:, 0:3]})
    ##store_data('jpar' + suffix, data={'x': timesb1, 'y': jparallel})
    ##store_data('jperp' + suffix, data={'x': timesb1, 'y': jperpvec[:, 0:3]})
    ##store_data('alpha' + suffix, data={'x': timesb1, 'y': alpha})
    ##store_data('alphaparallel' + suffix, data={'x': timesb1, 'y': alphaparallel})
    out_vars['timeseries' + suffix ] = timeseries
    out_vars['curlB' + suffix ] = curlB
    out_vars['divB' + suffix  ] = divB
    out_vars['jvec' + suffix] = jvec
    out_vars['jmag' + suffix] = jmag

    # Ignoring 'options' because they're specific to tplot an not used in this module rewrite.
    #
    ## set some options
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

    #return ['baryb', 'curlB', 'divB', 'jtotal', 'jpar', 'jperp', 'alpha', 'alphaparallel']
    return out_vars