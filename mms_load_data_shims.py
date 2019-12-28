"""
This module contains routines for loading MMS data


"""

from mms_curvature.mms_load_data import mms_load_data


# Leaving these imports listed here for my reference toward future implementation.
#from .fpi.mms_fpi_set_metadata import mms_fpi_set_metadata
#from .hpca.mms_hpca_set_metadata import mms_hpca_set_metadata
#from .feeps.mms_feeps_correct_energies import mms_feeps_correct_energies
#from .feeps.mms_feeps_flat_field_corrections import mms_feeps_flat_field_corrections
#from .feeps.mms_feeps_active_eyes import mms_feeps_active_eyes
#from .feeps.mms_feeps_split_integral_ch import mms_feeps_split_integral_ch
#from .feeps.mms_feeps_remove_bad_data import mms_feeps_remove_bad_data
#from .feeps.mms_feeps_remove_sun import mms_feeps_remove_sun
#from .feeps.mms_feeps_omni import mms_feeps_omni
#from .feeps.mms_feeps_spin_avg import mms_feeps_spin_avg

import re
import numpy as np 

def mms_load_fgm(trange=['2015-10-16', '2015-10-17'], probe='1', data_rate='srvy',
    level='l2', datatype='', varformat=None, prefix='', suffix='',
    keep_flagged=False, get_support_data=True, time_clip=False, no_update=False, available=False, notplot=False):
    """
    This function loads FGM data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for FGM include 'brst' 'fast' 'slow' 'srvy'. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            no datatype for FGM instrument (all science data are loaded)

        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a 
            "VAR_TYPE" attribute of "data".

        time_clip: bool
            Data will be clipped to the exact trange specified by the trange keyword.
            
        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        prefix: str
            The tplot variable names will be given this prefix.  By default, 
            no prefix is added.

        suffix: str
            The tplot variable names will be given this suffix.  By default, 
            no suffix is added.
            
    Returns:
        List of tplot variables created.

    """

    instrument='fgm'
    varformat_fetch = varformat

    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='fgm',
            datatype=datatype, varformat=varformat_fetch, prefix=prefix, suffix=suffix, get_support_data=get_support_data,
            time_clip=time_clip, no_update=no_update)
    
    #return tvars

    # the probes will need to be strings beyond this point
    if not isinstance(probe, list): probe = [probe]
    if not isinstance(data_rate, list): data_rate = [data_rate]
    if not isinstance(level, list): level = [level]
    
    probe = [('mms'+(str(p))) for p in probe]

    # remove flagged data
    if not keep_flagged:
        # Replacing this call by inlining the function contents.
        #mms_fgm_remove_flags(probe, data_rate, level, instrument, suffix=suffix)
        
        # From function's original comment string:
        #  Removes data flagged by the FGM 'flag' variable (flags > 0), 
        #  in order to only show science quality data by default.
        for this_probe in probe:
            for this_dr in data_rate:
                for this_lvl in level:
                    flag_var = this_probe+'_'+instrument+'_flag_'+this_dr+'_'+this_lvl+suffix
                    if flag_var in data.keys():
                        times, flags = data[flag_var].values()
                        flagged_data = np.where(flags != 0.0)[0]
            
                        for var_specifier in ['_b_gse_', '_b_gsm_', '_b_dmpa_', '_b_bcs_']:
                            var_name = this_probe+'_'+instrument+var_specifier+this_dr+'_'+this_lvl+suffix
                            if var_name in data.keys():
                                times, var_data = data[var_name].values()
                                var_data[flagged_data] = np.nan
        
        # Delete the flags variable if it was not originally requested
        if varformat is not None:
            regex = re.compile(varformat.replace("*", ".*"))
            datasets_to_delete = [dataset for dataset in data.keys() if not re.match(regex, dataset)]
            for dataset in datasets_to_delete:
                del data[dataset]
                if dataset in metadata.keys():
                    del metadata[dataset]

    return data,metadata



def mms_load_mec(trange=['2015-10-16', '2015-10-17'], probe='1', data_rate='srvy', 
    level='l2', datatype='ephts04d', varformat=None, prefix='', suffix='', get_support_data=False,
    time_clip=False, no_update=False, available=False, notplot=False):
    """
    This function loads MEC data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for MEC include ['brst', 'srvy']. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for MEC are: ['ephts04d', 'epht89q', 'epht89d']; default is 'ephts04d'

        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a 
            "VAR_TYPE" attribute of "data".

        time_clip: bool
            Data will be clipped to the exact trange specified by the trange keyword.
            
        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        prefix: str
            The tplot variable names will be given this prefix.  By default, 
            no prefix is added.

        suffix: str
            The tplot variable names will be given this suffix.  By default, 
            no suffix is added.

            
    Returns:
        List of tplot variables created.

    """

    data,metadata = mms_load_data(trange=trange, probe=probe, data_rate=data_rate, level=level, instrument='mec',
            descriptor=datatype, varformat=varformat, prefix=prefix, suffix=suffix, get_support_data=get_support_data,
            time_clip=time_clip, no_update=no_update)
    return data,metadata


