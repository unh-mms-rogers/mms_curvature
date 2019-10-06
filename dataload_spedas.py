import pyspedas
from pytplot import get_data
import logging

def DataLoad(trange=['2017-05-01', '2017-05-02/15:30:02'], data_rate='srvy', level='l2'):
    '''
    Loads all data needed for calculating magnnetic field curvature from MMS FGM data.  
    Uses pyspedas and pytplot.get_data for accessing the SDC API, file downloading, 
    data file version control, and CDF unpacking.

    Parameters:
    trange:     A list with two strings for the date range [tstart, tend]
                e.g. trange=['2017-05-01', '2017-05-02/15:30:02']

    data_rate:  The cadance of data which should be loaded.
                Options are 'srvy', 'brst'

    level:      The data level which will be loaded.  Use 'l2' unless you're sure otherwise.

    '''
    logging.info('Start DataLoad.')
    # load data files from SDC/local storage into tplot variables
    pyspedas.mms_load_mec(trange=trange, probe=['1', '2', '3', '4'], data_rate='srvy', level=level, time_clip=True)
    pyspedas.mms_load_fgm(trange=trange, probe=['1', '2', '3', '4'], data_rate=data_rate, level=level, time_clip=True)
    # extract data from tplot variables to numpy arrays.  NOTE: all done in GSM.
    postime1, pos1 = get_data('mms1_mec_r_gsm')
    postime2, pos2 = get_data('mms2_mec_r_gsm')
    postime3, pos3 = get_data('mms3_mec_r_gsm')
    postime4, pos4 = get_data('mms4_mec_r_gsm')
    magtime1, mag1 = get_data('mms1_fgm_b_gsm_'+data_rate+'_l2')
    magtime2, mag2 = get_data('mms2_fgm_b_gsm_'+data_rate+'_l2')
    magtime3, mag3 = get_data('mms3_fgm_b_gsm_'+data_rate+'_l2')
    magtime4, mag4 = get_data('mms4_fgm_b_gsm_'+data_rate+'_l2')
    # return all arrays
    logging.info('Returning from DataLoad.')
    return (postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4)

