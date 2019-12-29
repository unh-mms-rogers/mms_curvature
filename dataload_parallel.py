from . import mms_load_data_shims as shims
import logging

def DataLoad(trange=['2017-05-01', '2017-05-02/15:30:02'], data_rate='srvy', level='l2'):
    '''
    Loads all data needed for calculating magnnetic field curvature from MMS FGM data.  
    Uses a modified pymms for accessing the SDC API, file downloading, and file version control
    Uses a modified load_cdf from pyspedas for CDF unpacking.

    Parameters:
    trange:     A list with two strings for the date range [tstart, tend]
                e.g. trange=['2017-05-01', '2017-05-02/15:30:02']

    data_rate:  The cadance of data which should be loaded.
                Options are 'srvy', 'brst'

    level:      The data level which will be loaded.  Use 'l2' unless you're sure otherwise.

    '''
    
    logging.info('Start DataLoad2.')
    # load data files from SDC/local storage into tplot variables
    mec_data,mec_metadata = shims.mms_load_mec(trange=trange, probe=['1', '2', '3', '4'], data_rate='srvy', level=level, time_clip=True)
    fgm_data,fgm_metadata = shims.mms_load_fgm(trange=trange, probe=['1', '2', '3', '4'], data_rate=data_rate, level=level, time_clip=True)
    # extract data from tplot variables to numpy arrays.  NOTE: all done in GSM.
    postime1, pos1 = mec_data['mms1_mec_r_gsm'].values()
    postime2, pos2 = mec_data['mms2_mec_r_gsm'].values()
    postime3, pos3 = mec_data['mms3_mec_r_gsm'].values()
    postime4, pos4 = mec_data['mms4_mec_r_gsm'].values()
    magtime1, mag1 = fgm_data['mms1_fgm_b_gsm_'+data_rate+'_l2'].values()
    magtime2, mag2 = fgm_data['mms2_fgm_b_gsm_'+data_rate+'_l2'].values()
    magtime3, mag3 = fgm_data['mms3_fgm_b_gsm_'+data_rate+'_l2'].values()
    magtime4, mag4 = fgm_data['mms4_fgm_b_gsm_'+data_rate+'_l2'].values()
    # return all arrays
    logging.info('Returning from DataLoad2.')
    return (postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4)
