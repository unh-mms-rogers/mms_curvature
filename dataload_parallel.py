from . import mms_load_data_shims as shims
import logging

def DataLoad(trange=['2017-05-01', '2017-05-02/15:30:02'], data_rate='srvy', level='l2'):
    '''
    This function displays example usage of the data load methods found in mms_load_data_shims.py
    
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
    
    logging.info('Start parallel DataLoad.')
    # load data files from SDC/local storage into tplot variables
    mec_data,mec_metadata = shims.mms_load_mec(trange=trange, probe=['1', '2', '3', '4'], data_rate='srvy', level=level, time_clip=True)
    fgm_data,fgm_metadata = shims.mms_load_fgm(trange=trange, probe=['1', '2', '3', '4'], data_rate=data_rate, level=level, time_clip=True)
    logging.info('Returning from parallel DataLoad.')
    # This  is just a convenient structure for returning the imported data.
    return {
        "data": {
                    "mec": mec_data,
                    "fgm": fgm_data
                },
        "metadata": {
                    "mec": mec_metadata,
                    "fgm": fgm_metadata
                }
    }
