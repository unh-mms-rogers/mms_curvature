'''
calculate ion and electron gyroradius suing pyspedas

'''

import numpy as np
import mms_curvature as curvature

def DataLoadMoments(trange=['2017-05-01', '2017-05-02'], data_rate='srvy', level='l2', probe='1'):
    '''
    Loads all FPI moments data.

    Parameters:
    trange:     A list with two strings for the date range [tstart, tend]
                e.g. trange=['2017-05-01', '2017-05-02/15:30:02']

    data_rate:  The cadance of data which should be loaded.
                Options are 'srvy', 'brst'

    level:      The data level which will be loaded.  Use 'l2' unless you're sure otherwise.
    probe:      Define which probe data will be loaded from

    '''
    # Ensure data rate used is correct for FPI
    if data_rate == 'srvy': fpirate='fast'
    else: fpirate=data_rate

    # Load datafiles from SDC/local storage into pytplot variables
    fpidata,fpimeta = curvature.mms_load_fpi(trange=trange, probe=probe, data_rate=fpirate, level=level, datatype=['dis-moms', 'des-moms'], time_clip=True)
    
    distime, distempperp = fpidata['mms'+probe+'_dis_tempperp_'+fpirate].values()
    destime, destempperp = fpidata['mms'+probe+'_des_tempperp_'+fpirate].values()
    #extract data from tplot variables
    #distime, distempperp = get_data('mms'+probe+'_dis_tempperp_'+fpirate)
    #destime, destempperp = get_data('mms'+probe+'_des_tempperp_'+fpirate)
    
    # return all arrays
    return (distime, distempperp, destime, destempperp)

def CalcRadius(part_time=None, part_tempperp=None, b_time=None, b_mag=None, part_mass=1.672622e-27, part_q=1.602177e-19):
    '''
    Calculates gyroradius of a given particle.  Interpolates particle data
    to magnetic field cadance.

    Parameters:
    part_time:      Time series of the particle moment data for interpolation
    part_tempperp:  Perpendicular temperature of particle aligned with part_time
    b_time:         Time series of the magnetic field data
    b_mag:          Magnitude of the magnetic field aligned with b_time
    part_mass:      Mass of the particle to eb calculated in kg 
                    (defualts to proton mass)
    part_q:         Magnitude of charge of the particle to be calculated in Coulombs 
                    (defaults to elementary charge)
    '''
    
    # Interpolate particle data to megnetic field cadance and time series
    tperp = np.interp(b_time, part_time, part_tempperp)
    
    # define scaling constants.  Assumes tperp is in eV, bmag is in nT
    constants = np.sqrt(2*part_mass*11605*1.38649e-23)/(1e-9*part_q)
    
    # Calculate gyroradius
    r_g = constants * np.sqrt(tperp)/b_mag
    
    # Return gyroradius in meters
    return (r_g)


