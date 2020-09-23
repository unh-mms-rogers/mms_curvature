'''
python3 + numpy routine to calculate magnetic curvature from MMS data.

Uses pyspedas for MMS data file loading

*****************************************************************************
Example script for loading required data and calculating the curvature vector:
    
import time
import re
import numpy as np
from mms_curvature import DataLoad, mms_Grad, mms_Curvature

timeStart =  time.strftime("%H:%M:%S", time.localtime())
print("Files Loading:")

data,metadata = DataLoad(trange=['2017-05-04', '2017-05-05'])
postimes = [data['mec'][n]['x'] for n in data['mec'] if re.compile('mms\d_mec_r_gsm_srvy_l2').match(n)]
posvalues = [data['mec'][n]['y'] for n in data['mec'] if re.compile('mms\d_mec_r_gsm_srvy_l2').match(n)]
magtimes = [data['fgm'][n]['x'] for n in data['fgm'] if re.compile('mms\d_fgm_b_gsm_srvy_l2').match(n)]
magvalues = [data['fgm'][n]['y'] for n in data['fgm'] if re.compile('mms\d_fgm_b_gsm_srvy_l2').match(n)]

print("Time started: ", timeStart)
print("Time Loaded: ", time.strftime("%H:%M:%S", time.localtime()))

grad_Harvey, bm, bmag, rm, t_master = mms_Grad(postimes, posvalues, magtimes, magvalues)
curve_Harvey = mms_Curvature(grad_Harvey, bm)

np.savetxt("t_master.csv", t_master, delimiter=",")
np.savetxt("curve_Harvey.csv", curve_Harvey, delimiter=",")
np.save("grad_Harvey.npy", grad_Harvey)

# end
****************************************************************************

---TJR 10.09.2020


NOTE: Need to update the docstring for mms_Grad after refactoring
'''
import numpy as np

def mms_Grad(postimes=None, posvalues=None, magtimes=None, magvalues=None, normalize=True):
    
    '''
    Calculates spacial gradient and curvature vector of the magnetic field.  
    Returns those and the master time (interpolated from FGM data) as numpy arrays

    normalize:  if True normalizes magnetic field vectors before continusing
                calculation; required for calculating curvature

                if False leaves magnetic field as full vector with magnitude;
                required for calculating curl and divergence
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
        for bird in range(numBirds):
            bn[bird] = magvalues[bird][:,0:3]/magvalues[bird][:,3,np.newaxis]
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
    for bird in range(numBirds):
        Bvals[bird] = np.interp(t_master, magtimes[bird], magvalues[bird][:,3])
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
    
    # calculate position and magnetic field at mesocenter of the fleet
    
    
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
    
    # For each timestep t, dtemp[:,t,:] now looks like this:
    # dtemp[:,t] = [[dB_x*dR_x, dB_y*dR_y, dB_z*dR_z],
    #               [dB_x*dR_y, dB_y*dR_z, dB_z*dR_x],
    #               [dB_x*dR_z, dB_y*dR_x, dB_z*dR_y]]
    
    # The below constructs dbdr by twisting the dtemp array to properly place the diagonals for dbdr.
    #
    # dbdr[t] = [[dtemp[0,t,0], dtemp[1,t,0], dtemp[2,t,0],
    #            [dtemp[2,t,1], dtemp[0,t,1], dtemp[1,t,1],
    #            [dtemp[1,t,2], dtemp[2,t,2], dtemp[0,t,2]]
    #
    # ===
    #
    # dbdr[t] = [[dB_x*dR_x, dB_x*dR_y, dB_x*dR_z],
    #            [dB_y*dR_x, dB_y*dR_y, dB_y*dR_z],
    #            [dB_z*dR_x, dB_z*dR_y, dB_z*dR_z]]
    #
    
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

    tmpHarv = np.ndarray((3, t_master.shape[0], 3))
    tmpHarv[0] = np.divide(np.einsum('...i,...i',np.moveaxis(Rinv,1,-1),dbdr),np.square(numBirds))
    tmpHarv[1] = np.divide(np.einsum('...i,...i',np.moveaxis(Rinv,1,-1),np.roll(dbdr,-1,1)),np.square(numBirds))
    tmpHarv[2] = np.divide(np.einsum('...i,...i',np.moveaxis(Rinv,1,-1),np.roll(dbdr,-2,1)),np.square(numBirds))
    
    # Constructs the gradient matrix for each timestep from the diagonals calculated in the above steps.
    grad_Harvey = \
                  np.einsum('...i,...ij->...ij',tmpHarv[0],np.identity(3)) + \
                  np.einsum('...i,...ij->...ij',tmpHarv[1],(np.roll(np.identity(3),-1,axis=0))) + \
                  np.einsum('...i,...ij->...ij',tmpHarv[2],(np.roll(np.identity(3),-2,axis=0)))

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
    

    return grad_Harvey, bm, bmag, rm, t_master




def mms_Curvature(grad, bm):
    '''
    function to calculate magnetic field line curvature vector k = b · grad(b) from the
    magnetic field spacial gradient (grad) and the normalized magnetic field vector at 
    the barycenter of the spacecraft formation (bm)

    Inputs:
    grad    : a time series array with dimensions t x 3 x 3 representing the spacial
              gradient of the normalized magnetic field grad(b) at each time step.
              Assumed to be the output from the mms_Grad function described in this 
              library module.

    bm      : a time series array with dimensions t x 3 representing the normalized 
              vector magnetic field b = B/|B| at each time step.  Assumed to be the
              output of the mms_Grad function described in this library.

    ---NB: 'grad' and 'bm' are assumed to have identical timesteps, beginning, and 
           ending times, as expected from the outputs of the mms_Grad function 
           described in this library.

    Outputs:
    curve_Harvey    : a time series array with dimensions t x 3 representing the 
                      magnetic field line curvature vector in 3 dimensions at the
                      same time steps as used for the input arrays.
    '''


    # And now the final curvature may be calculated by simple matrix multiplication for each timestep.
    # The below gives identical results as, but is much more efficient than: 
    #   for t in range(t_master.shape[0]): curve_Harvey[t] = np.matmul(grad_Harvey[t], bm[t])
    curve_Harvey = np.einsum('...ij,...j', grad, bm)
    
    return curve_Harvey

    

def mms_CurlB(Grad):
    '''
    Function to calculate the curl of the magnetic field by applying a rank 3 
    Levi-Civita tensor to the spacial gradient of the full vector magnetic 
    field (NOTE: NOT the gradient of the normalized vector magnetic field as
    used in the curvature vector calculation).

    Inputs:
    Grad    : A time series array with dimensions t x 3 x 3 representing the 
              spacial gradient of the vector magnetic field grad(B) at each
              time step.  Assumed to be the output from the mms_Grad function
              described in this library module.

    Outputs:
    CurlB   : A time series array with dimensions t x 3 representing the curl
              of the vector magnetic field in 3 dimensions at the same time
              steps as used for hte input Grad array.
    '''
    
    # Define the rank 3 Levi-Civita tensor
    LevCiv3 = np.zeros((3,3,3))
    LevCiv3[0,1,2] = LevCiv3[1,2,0] = LevCiv3[2,0,1] = 1
    LevCiv3[0,2,1] = LevCiv3[2,1,0] = LevCiv3[1,0,2] = -1

    CurlB = np.einsum('ijk,...jk',LevCiv3, Grad)

    return CurlB



def mms_DivB(Grad):
    '''
    Function to calculate the divergence of the magnetic field by taking the
    trace of the spacial gradient fo the full vector magnetic field (NOTE: 
    NOT the gradient of the normalized vector magnetic foeld as used in the 
    curvature vector calculation).

    Inputs:
    Grad    : A time series array with dimensions t x 3 x 3 representing the 
              spacial gradient of the vector magnetic field grad(B) at each
              time step.  Assumed to be the output from the mms_Grad function
              described in this library module.

    Outputs:
    DivB    : A time series array with dimension t representing the divergence
              of the vector magnetic field div(B) at the same time steps as 
              used for the input Grad array.
    '''

    DivB = np.einsum('...ii', Grad)

    return DivB
    





