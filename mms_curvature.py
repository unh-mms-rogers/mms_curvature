'''
python3 + numpy routine to calculate magnetic curvature from MMS data.

Uses pyspedas for MMS data file loading

*****************************************************************************
Example script for loading required data and calculating the curvature vector:
    
import time
import numpy as np
from mms_curvature import DataLoad, Curvature

timeStart =  time.strftime("%H:%M:%S", time.localtime())
print("Files Loading:")

postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4 = DataLoad(trange=['2017-05-04', '2017-05-05'])

print("Time started: ", timeStart)
print("Time Loaded: ", time.strftime("%H:%M:%S", time.localtime()))

t_master, grad_Harvey, curve_Harvey = Curvature(postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4)

np.savetxt("t_master.csv", t_master, delimiter=",")
np.savetxt("curve_Harvey.csv", curve_Harvey, delimiter=",")
np.save("grad_Harvey.npy", grad_Harvey)

# end
****************************************************************************

---AJR 08.22.2019
'''
import numpy as np

def Curvature(postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4, report_all=False, with_uncertainty=False): 
    '''
    Calculates spacial gradient and curvature vector of the magnetic field.  
    Returns those and the master time (interpolated from FGM data) as numpy arrays
    '''
    
    # normalize magnetic fields
    bn1 = mag1[:,0:3]/mag1[:,3,np.newaxis]
    bn2 = mag2[:,0:3]/mag2[:,3,np.newaxis]
    bn3 = mag3[:,0:3]/mag3[:,3,np.newaxis]
    bn4 = mag4[:,0:3]/mag4[:,3,np.newaxis]

    # find probe with latest beginning point for magnetic field data
    mastersc = np.argmax([magtime1[0], magtime2[0], magtime3[0], magtime4[0]])
    # find earliest end time for all space craft in range
    tend = min(magtime1[-1], magtime2[-1], magtime3[-1], magtime4[-1])
    tend_i = 0  # initialize counting index for finding ending index for master S/C
    
    # initialize master time sequence by trimming time from last S/C to start (mastersc)
    if mastersc == 0:
        while magtime1[tend_i] < tend: tend_i += 1
        t_master = magtime1[0:(tend_i+1)]
    elif mastersc == 1:
        while magtime2[tend_i] < tend: tend_i += 1
        t_master = magtime2[0:(tend_i+1)]
    elif mastersc == 2:
        while magtime3[tend_i] < tend: tend_i += 1
        t_master = magtime3[0:(tend_i+1)]
    elif mastersc == 3:
        while magtime4[tend_i] < tend: tend_i += 1
        t_master = magtime4[0:(tend_i+1)]
    
    # master mag field data arr, with interpolated values
    # Magnetic field data, interpolated to the previously determined master time sequence
    barr=np.ndarray((4,t_master.shape[0],3))
    
    barr[0,:,0] = np.interp(t_master, magtime1, bn1[:,0]) # MMS1
    barr[0,:,1] = np.interp(t_master, magtime1, bn1[:,1]) 
    barr[0,:,2] = np.interp(t_master, magtime1, bn1[:,2]) 
    
    barr[1,:,0] = np.interp(t_master, magtime2, bn2[:,0]) # MMS2
    barr[1,:,1] = np.interp(t_master, magtime2, bn2[:,1]) 
    barr[1,:,2] = np.interp(t_master, magtime2, bn2[:,2]) 
    
    barr[2,:,0] = np.interp(t_master, magtime3, bn3[:,0]) # MMS3
    barr[2,:,1] = np.interp(t_master, magtime3, bn3[:,1]) 
    barr[2,:,2] = np.interp(t_master, magtime3, bn3[:,2]) 
    
    barr[3,:,0] = np.interp(t_master, magtime4, bn4[:,0]) # MMS4
    barr[3,:,1] = np.interp(t_master, magtime4, bn4[:,1]) 
    barr[3,:,2] = np.interp(t_master, magtime4, bn4[:,2]) 


    # Calculate average |B| for export
    B1 = np.interp(t_master, magtime1, mag1[:,3])
    B2 = np.interp(t_master, magtime2, mag2[:,3])
    B3 = np.interp(t_master, magtime3, mag3[:,3])
    B4 = np.interp(t_master, magtime4, mag4[:,3])
    bmag = np.average([B1, B2, B3, B4], axis=0)



    # The below code is preseved in case a future need calls for timesteps to by trimmed
    #  in this function instead of during data load of input for this function.
    '''
    # need to clear extranious data points from positional data for
    #   np.interp to work as intended.
    bcnt = 0        # MMS1
    while postime1[bcnt] < t_master[0] and bcnt: bcnt = bcnt + 1
    ecnt = bcnt
    while postime1[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos1 = pos1[(bcnt-1):(ecnt+1),:]
    postime1 = postime1[(bcnt-1):(ecnt+1),:]
   
    bcnt = 0        # MMS2
    while postime2[bcnt] < t_master[0]: bcnt = bcnt + 1
    ecnt = bcnt
    while postime2[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos2 = pos2[(bcnt-1):(ecnt+1),:]
    postime2 = postime2[(bcnt-1):(ecnt+1),:]

    bcnt = 0        # MMS3
    while postime3[bcnt] < t_master[0]: bcnt = bcnt + 1
    ecnt = bcnt
    while postime3[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos3 = pos3[(bcnt-1):(ecnt+1),:]
    postime3 = postime3[(bcnt-1):(ecnt+1),:]
    
    bcnt = 0        # MMS4
    while postime4[bcnt] < t_master[0]: bcnt = bcnt + 1
    ecnt = bcnt
    while postime4[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos4 = pos4[(bcnt-1):(ecnt+1),:]
    postime4 = postime4[(bcnt-1):(ecnt+1),:]
    ''' 
    
    # master position data array, with interpolated value
    # Spacecraft position data, interpolated to the previously determined master time sequence
    rarr = np.ndarray((4,t_master.shape[0],3))
    
    rarr[0,:,0] = np.interp(t_master, postime1, pos1[:,0]) # MMS1
    rarr[0,:,1] = np.interp(t_master, postime1, pos1[:,1])
    rarr[0,:,2] = np.interp(t_master, postime1, pos1[:,2])
    
    rarr[1,:,0] = np.interp(t_master, postime2, pos2[:,0]) # MMS2
    rarr[1,:,1] = np.interp(t_master, postime2, pos2[:,1])
    rarr[1,:,2] = np.interp(t_master, postime2, pos2[:,2])
    
    rarr[2,:,0] = np.interp(t_master, postime3, pos3[:,0]) # MMS3
    rarr[2,:,1] = np.interp(t_master, postime3, pos3[:,1])
    rarr[2,:,2] = np.interp(t_master, postime3, pos3[:,2])
    
    rarr[3,:,0] = np.interp(t_master, postime4, pos4[:,0]) # MMS4
    rarr[3,:,1] = np.interp(t_master, postime4, pos4[:,1])
    rarr[3,:,2] = np.interp(t_master, postime4, pos4[:,2])
    
    # Now all magnetic fields and positional data of of the same cadence and at the same times for each index
    # Indices are: [s/c(0=mms1, 1=mms2, 2=mms3, 3=mms4), time_step, vector(0=x, 1=y, 2=z)]
    # ie.  rarr[<spacecraft>, <timestep_index>, <cartesian_component_of_vector>]
    # eg.  Y-position of mms4 at first time step:  rarr[3,0,1]
    
    # calculate position and normalized magnetic field at mesocenter of the fleet
    
    
    # Commenting old implementation
    #rm = np.ndarray((t_master.shape[0], 3))
    #for i in range(t_master.shape[0]):      # populate each element of the mesocenter with the average across all four s/c
    #    for j in range(3):                  # there's got to be a way to vectorize this
    #        #print("rm:", rm.shape, "rarr:", rarr.shape)
    #        rm[i,j] = (1./4.)*(rarr[0,i,j]+rarr[1,i,j]+rarr[2,i,j]+rarr[3,i,j])
    #
    #bm = np.ndarray((t_master.shape[0], 3))
    #for i in range(t_master.shape[0]):      # populate each element of the mesocenter with the average across all four s/c
    #    for j in range(3):                  # there's got to be a way to vectorize this
    #        bm[i,j] = (1./4.)*(barr[0,i,j]+barr[1,i,j]+barr[2,i,j]+barr[3,i,j])
    # End of old implementation.
    
    # Vectorized version of above, using numpy built-in ufuncs.
    # - inner operation `np.add.reduce(input, axis=0)`: .reduce collapses the input array along the specified
    #   axis (0 by default), using the prefaced ufunc (add) to merge array values.
    #     Layman's version:  Add up the respective coordinate components from each bird.
    # - outer operation `np.divide(input1,input2)`: Divides each element in input1 by the respective element
    #   in input 2.  If input2 is of lesser dimensions than input1, input2 will be expanded/broadcast to fit.
    #     Layman's version:  Divide each component of the inner results by the number of birds.
    rm = np.divide(np.add.reduce(rarr),rarr.shape[0])
    bm = np.divide(np.add.reduce(barr),barr.shape[0])
    
    # Calculate volumetric tensor (Harvey, Ch 12.4, Eq 12.23, from "Analysis Methods for Multi-Spacecraft Data" Paschmann ed.)
    
    #Rvol = np.ndarray((t_master.shape[0], 3, 3))
    #for i in range(t_master.shape[0]):
    #    #rvol_step = np.zeros([3,3])    # Stepwise method gives same result as explicit below
    #    #for sc in range(4):
    #    #    rvol_step = rvol_step + np.outer(rarr[sc,i,:], rarr[sc,i,:])
    #    ## endfor
    #    #Rvol[i,:,:] = (1./4.) * (rvol_step - np.outer(rm[i,:], rm[i,:]))
    #    #Rvol[i,:,:] = (1./4.)*((np.outer(rarr[0,i,:], rarr[0,i,:]) + np.outer(rarr[1,i,:], rarr[1,i,:]) + np.outer(rarr[2,i,:], rarr[2,i,:]) + np.outer(rarr[3,i,:], rarr[3,i,:])) - np.outer(rm[i,:], rm[i,:]))     # give same result as stepwise above
    
    # Intermediate variable to hold the self-outer-product of rm at each timestep
    rmOuter = np.einsum('...i,...j->...ij',rm,rm)  # explicit form 'outer product' use of EinSum, broadcast across leading dimensions
    
    # Intermediate variable to hold the self-outer-product of rarr, per bird, at each timestep
    rarrOuter = np.einsum('...i,...j->...ij',rarr,rarr)  # Same as line above
    
    Rvol = np.divide(np.add.reduce(rarrOuter) - rmOuter, rarr.shape[0])     # give same result as stepwise Rvol code commented out above.
    # Brief description of operations in above line:
    #  All of the following occur for each timestep...
    #  1)  Collapse the self-outer-products of rarr by summing across all spacecraft
    #  2)  From above result, subtract the self-outer-product of rm (position of mesocenter)
    #  3)  Divide resultant array from above step by the number of spacecraft
    
    # Pre-calculate the inverse array of Rvol here to avoid needless recalculation later.
    Rinv = np.linalg.inv(Rvol)
    
    # Stepwise calculation of gradient and curvature using the Harvey method
    #a_ne_b_list=[[1,2,3],[2,3],[3]]
    #grad_Harvey = np.ndarray((t_master.shape[0], 3, 3))
    #curve_Harvey = np.ndarray((t_master.shape[0], 3))    
    #for t in range(t_master.shape[0]):      # steps through each time step
    #    for i in range(3):                  # for final i-component of the gradient 
    #        for j in range(3):              # for final j-component of the gradient
    #            dbdr = np.zeros(3)          # a != b summation row vector from Harvey.  Re-initialized as zeros for each i,j
    #            for k in range(3):          # step through k-index.  May be able to eliminate this with vectorization later
    #                for a in range(3):      # step through spacecraft MMS1-3; MMS4 done implicitly
    #                    for b in a_ne_b_list[a]:    # Does not contain MMS1 (done by stepping through it in previous loop); provides sc_a != sc_b summation in Harvey (12.18)
    #                        dbdr[k] = dbdr[k] + ((barr[a,t,i] - barr[b,t,i]) * (rarr[a,t,k] - rarr[b,t,k]))
    #                    # endfor
    #                # endfor
    #            # endfor
    #            # grad_Harvey[t,i,j] = (1./16.) * np.matmul(dbdr, Rinv[t,:,j])     # Gives the same result as below
    #            grad_Harvey[t,i,j] = (1./16.) * np.matmul(Rinv[t,:,j], dbdr)     # Gives the same result as below
    #            #grad_Harvey[t,i,j] = (1./16.) * np.matmul(dbdr, np.linalg.inv(Rvol[t,:,:])[:,j])    # Maybe linalg.inv doesn't vectorize the way I think?
    #        # endfor
    #    # curve_Harvey[t,:] = np.matmul(bm[t,:], grad_Harvey[t,:,:])               # same thing, probably just should be matrix mult.
    #    curve_Harvey[t,:] = np.matmul(grad_Harvey[t,:,:], bm[t,:])               # Order of matmul has BIG effect!
    ## endfor
    
    # Vectorized matrix operations to calculate the above.  Saves a lot of compute time at the expense of a little memory.
    tmpR = np.repeat(rarr[np.newaxis,:,:,:],rarr.shape[0],axis=0)  # Stretch the array to be 2-D instead of 1-D for the sats.  Required for next operation.
    triR = np.triu(np.rollaxis(np.rollaxis(np.transpose(tmpR,axes=(1,0,2,3)) - tmpR, -1), -1))  # This produces a triangular matrix of dR=D_a - D_b, for all [a != b]
    
    tmpB = np.repeat(barr[np.newaxis,:,:,:],barr.shape[0],axis=0)  # Same as above, but with B instead
    triB = np.triu(np.rollaxis(np.rollaxis(np.transpose(tmpB,axes=(1,0,2,3)) - tmpB, -1), -1)) # Again, dB=B_a - B_b, for all [a != b]
    
    #Example of effect of above operations:
    # Each b_i below is a 3-vector
    # 
    # Line 249 (at each timestep):
    # tmpB = 
    # [[b_1, b_2, b_3, b_4],
    #  [b_1, b_2, b_3, b_4],
    #  [b_1, b_2, b_3, b_4],
    #  [b_1, b_2, b_3, b_4]]
    # 
    # Line 250, two steps, first array is intermediate form (again, at each timestep):
    # [[b_1-b_1,  b_1-b_2,  b_1-b_3, b_1-b_4],
    #  [b_2-b_1,  b_2-b_2,  b_2-b_3, b_2-b_4],
    #  [b_3-b_1,  b_3-b_2,  b_3-b_3, b_3-b_4],
    #  [b_4-b_1,  b_4-b_2,  b_4-b_3, b_4-b_4]]
    # 
    # triB =
    # [[0      ,  b_1-b_2,  b_1-b_3, b_1-b_4],
    #  [0      ,  0      ,  b_2-b_3, b_2-b_4],
    #  [0      ,  0      ,  0      , b_3-b_4],
    #  [0      ,  0      ,  0      , 0      ]]
    
    # Calculate the partial components for dbdr
    dtemp = np.ndarray((3, t_master.shape[0], 3))
    dtemp[0] = np.einsum('...ab,...ab',triB,triR) #This gets us the diagonals of dbdr for B_i and R_i (eg, both x components, both y, ...)
    dtemp[1] = np.einsum('...ab,...ab',triB,np.roll(triR,-1,axis=1)) #This gets us the diagonals of dbdr for B_i and R_mod(i+1)  (eg, B_x * R_y, ...)
    dtemp[2] = np.einsum('...ab,...ab',triB,np.roll(triR,-2,axis=1)) #This gets us the diagonals of dbdr for B_i and R_mod(i+2)  (eg, B_y * R_x, ...)
    
    # Constructs dbdr matrix for each timestep, where dbdr[i,j] will return the relavant dB_i*dR_j resultant vector
    dbdr = np.einsum('...i,...ij->...ij',dtemp[0],np.identity(3)) + \
            np.einsum('...i,...ij->...ij',dtemp[1],(np.roll(np.identity(3),-1,axis=0))) + \
            np.einsum('...i,...ij->...ij',dtemp[2],(np.roll(np.identity(3),-2,axis=0)))
    
    # This calculates and holds the diagonals for the Harvey gradient.  I'm sure there's some simpler way to calculate this, but I haven't found it yet.
    # This eventually gets us to the Harvey gradients in the same manner as we got dbdr above.

    numBirds = 4    # Set for MMS.  Can be changed with rest of code to expand for a larger fleet

    tmpHarv = np.ndarray((3, t_master.shape[0], 3))
    tmpHarv[0] = np.divide(np.einsum('...i,...i',np.moveaxis(Rinv,1,-1),dbdr),np.square(numBirds))
    tmpHarv[1] = np.divide(np.einsum('...i,...i',np.moveaxis(Rinv,1,-1),np.roll(dbdr,-1,1)),np.square(numBirds))
    tmpHarv[2] = np.divide(np.einsum('...i,...i',np.moveaxis(Rinv,1,-1),np.roll(dbdr,-2,1)),np.square(numBirds))
    
    # Constructs the gradient matrix for each timestep from the diagonals calculated in the above steps.
    grad_Harvey = np.transpose( \
                  np.einsum('...i,...ij->...ij',tmpHarv[0],np.identity(3)) + \
                  np.einsum('...i,...ij->...ij',tmpHarv[1],(np.roll(np.identity(3),-1,axis=0))) + \
                  np.einsum('...i,...ij->...ij',tmpHarv[2],(np.roll(np.identity(3),-2,axis=0))) \
                  , (0,2,1)) # Due to all the matrix rolling shenanigans leading up to this, we need to transpose the matrix from each timestep.

    # And now the final curvature may be calculated by simple matrix multiplication for each timestep.
    # The below gives identical results as, but is much more efficient than: 
    #   for t in range(t_master.shape[0]): curve_Harvey[t] = np.matmul(grad_Harvey[t], bm[t])
    curve_Harvey = np.einsum('...ij,...j', grad_Harvey, bm)
    
    ## List of references for how numpy.einsum operates:
    # 
    # Official docs on einsum (as of numpy 1.18): https://numpy.org/doc/1.18/reference/generated/numpy.einsum.html
    # General explanation of einsum: https://stackoverflow.com/a/33641428
    # Levi-Cevita and einsum: https://stackoverflow.com/a/20890335
    # A specific instance of Levi-Cevita with einsum:  https://stackoverflow.com/a/20910319
    
    # Solenoid correction was implented then removed due to negligable impact on results.
    # Original stepwise code for solenoid correction is preseved in the below string in case of future need.
    '''
    # Solenoid correction from Harvey (12.20)
    # lm = np.ndarray((t_master.shape[0]))
    lm = np.divide(np.trace(grad_Harvey, axis1=1, axis2=2), np.trace(Rinv, axis1=1, axis2=2))
    
    nRinv = np.ndarray(Rinv.shape)
    for t in range(t_master.shape[0]):
        nRinv[t,:,:] = lm[t]*Rinv[t,:,:]

    sol_grad_Harvey = grad_Harvey - nRinv
    sol_curve_Harvey = np.ndarray((t_master.shape[0], 3))
    for t in range(t_master.shape[0]):
        sol_curve_Harvey[t,:] = np.matmul(sol_grad_Harvey[t,:,:], bm[t,:])
    '''         # Solenoid correction has little effect, which is not surprising

    outputs = (t_master, grad_Harvey, curve_Harvey)

    if with_uncertainty:
        float_max = np.finfo(float).max
        '''
        Breakdown of the following pair of lines:
        np.minimum.reduce(...) -    This will be producing a result array of rank 2 (timesteps, spacial direction 3-vector)
                                    from a rank 4 input array, by returning the minimum value from the num_birds X num_birds
                                    sub-matrix containing value differences between each pair of spacecraft.
        np.abs(input_array) -   We only care about relative distance of the values between spacecraft, so we compare the
                                absolute values of <input_array>.
        axis=(-1,-2) -  Input for the reduce function.  Collapse the last (-1) and next-to-last (-2) dimensions of the input.
        initial=float_max  -  Sets the initial value to compare with for each set to max possible float value.
                              This ensures that any comparison will result in the recorded value.
        where=np.logical_not(np.isclose(triR,0,atol=0)) - 
                        This means we're only performing comparisons with non-zero values.  In the vanishingly rare
                        circumstance that the smallest difference between two spacecraft _is_ 0, we use the
                        second-smallest difference.
                        If there is somehow a measurement where all spacecraft perfectly agree, the min value will
                        be set to 0 by the two lines following the initial reduction.
        '''
        minR = np.minimum.reduce(np.abs(triR), axis=(-1,-2), initial=np.finfo(float).max, where=np.logical_not(np.isclose(triR,0,atol=0)))
        minB = np.minimum.reduce(np.abs(triB), axis=(-1,-2), initial=np.finfo(float).max, where=np.logical_not(np.isclose(triB,0,atol=0)))
        minR[minR==float_max] = 0 # Resets any timesteps for which all differences were 0.  (Ludicriously unlikely, but should still do the check.)
        minB[minB==float_max] = 0
        outputs += (minR, minB)

        ###  Inserting lines for calculating uncertainty of the gradient tensor

        from mms_curvature.utils.uncertainty import sigmaAB

        uB = 0.1    # default instrument uncertainty of the magnetometer is 0.1 nT
        uR = 0.1    # default measurement uncertainty for position is 100m = 0.1km (I'd love to make this better)

        # ugrad = array of uncertainty ratios (sigma_f/f) for each element of the grad_b array, assuming the smallest separation in both B and R and so the dominant source of error, maximum estimate

        ugrad=np.array([[sigmaAB(minB[:,0], uB, minR[:,0], uR, op='*')[0], sigmaAB(minB[:,0], uB, minR[:,1], uR, op='*')[0], sigmaAB(minB[:,0], uB, minR[:,2], uR, op='*')[0]] \ 
        ,[sigmaAB(minB[:,1], uB, minR[:,0], uR, op='*')[0], sigmaAB(minB[:,1], uB, minR[:,1], uR, op='*')[0], sigmaAB(minB[:,1], uB, minR[:,2], uR, op='*')[0]] \ 
        ,[sigmaAB(minB[:,2], uB, minR[:,0], uR, op='*')[0], sigmaAB(minB[:,2], uB, minR[:,1], uR, op='*')[0], sigmaAB(minB[:,2], uB, minR[:,2], uR, op='*')[0]]])

        # umag = uncertainty ratio (sigma_f/f) for each vector component of the magnetic field at the mesocenter

        umag = np.divide(uB, np.multiply(bm, bmag))

        # uk = uncertainty ratio (sigma_f/f) for the curvature vector 'k'

        uk = np.einsum('...ij,...i', np.square(ugrad), np.square(umag))

        outputs += uk





        # first input assumed to be magnetic field vector; second assumed to be position



    if report_all:
        # Use barr and B* to find minimum dB along each axis
        a_ne_b_list=[[1,2,3],[2,3],[3]]
        B_arr = [B1, B2, B3, B4]
        dBmin = np.full((t_master.shape[0],3), np.inf)     # minimum dB for x,y,z GSM 
        for t in range(t_master.shape[0]):          # step through the time series
            for a in range(3):                      # a = MMS 1(0)  -> 3(2); MMS 4 done implicitly
                for b in a_ne_b_list[a]:            # b for sc_a != sc_b to do all possible combonations
                    tmp_x = (barr[a,t,0] * B_arr[a][t]) - (barr[b,t,0] * B_arr[b][t])
                    tmp_y = (barr[a,t,1] * B_arr[a][t]) - (barr[b,t,1] * B_arr[b][t])
                    tmp_z = (barr[a,t,2] * B_arr[a][t]) - (barr[b,t,2] * B_arr[b][t])
                    if tmp_x < dBmin[t,0]: dBmin[t,0] = tmp_x
                    if tmp_y < dBmin[t,1]: dBmin[t,1] = tmp_y   # finds minimum change in B along each axis
                    if tmp_z < dBmin[t,2]: dBmin[t,2] = tmp_z
                #end for
            #end for
        #end for

        outputs += (rarr, barr, rm, bm, bmag, dBmin, Rvol, Rinv)  #used for troubleshooting

    return outputs
    





# end
