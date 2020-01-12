"""
This module contains routines for loading MMS data


"""

from .mms_load_data import mms_load_data


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

def mms_load_fgm(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy',
        level='l2',
        datatype='',
        varformat=None,
        prefix='',
        suffix='',
        keep_flagged=False,
        get_support_data=True,
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
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
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to True.

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

    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument=instrument,
            descriptor=datatype, varformat=varformat, prefix=prefix, suffix=suffix, get_support_data=get_support_data,
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



def mms_load_mec(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy', 
        level='l2',
        datatype='ephts04d',
        varformat=None,
        prefix='',
        suffix='',
        get_support_data=False,
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
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
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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



def mms_load_hpca(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy', 
        level='l2',
        datatype='moments',
        get_support_data=True,
        time_clip=False,
        no_update=False,
        varformat=None,
        prefix='',
        suffix='',
        center_measurement=False,
        available=False,
        notplot=False ):
    """
    This function loads HPCA data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for HPCA include 'brst', 'srvy'. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for HPCA are 'moments' and 'ion'; the default is 'moments'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to True.

        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        time_clip: bool
            Data will be clipped to the exact trange specified by the trange keyword.
            
        prefix: str
            The tplot variable names will be given this prefix.  By default, 
            no prefix is added.

        suffix: str
            The tplot variable names will be given this suffix.  By default, 
            no suffix is added.

        center_measurement: bool
            If True, the CDF epoch variables are time-shifted to the middle
            of the accumulation interval by their DELTA_PLUS_VAR and
            DELTA_MINUS_VAR variable attributes

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)

            
    Returns:
        List of tplot variables created.

    """

    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='hpca',
            descriptor=datatype, varformat=varformat, prefix=prefix, suffix=suffix, get_support_data=get_support_data,
            time_clip=time_clip, no_update=no_update, center_measurement=center_measurement)
    
    return data,metadata



def mms_load_fpi(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='fast',
        level='l2',
        datatype=['des-moms', 'dis-moms'],
        varformat=None,
        prefix='',
        suffix='',
        get_support_data=False,
        time_clip=False,
        no_update=False,
        center_measurement=False,
        available=False,
        notplot=False ):
    """
    This function loads FPI data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for FPI include 'brst', 'fast'. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for FPI are:
             'des-moms', 'dis-moms' (default)
             'des-dist', 'dis-dist'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        center_measurement: bool
            If True, the CDF epoch variables are time-shifted to the middle
            of the accumulation interval by their DELTA_PLUS_VAR and
            DELTA_MINUS_VAR variable attributes

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)


            
    Returns:
        List of tplot variables created.

    """

    data,metadata = mms_load_data(trange=trange, probe=probe, data_rate=data_rate, level=level, instrument='fpi',
            descriptor=datatype, varformat=varformat, prefix=prefix, suffix=suffix, get_support_data=get_support_data,
            time_clip=time_clip, no_update=no_update, center_measurement=center_measurement, notplot=notplot)
    
    return data,metadata



def mms_load_scm(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy', 
        level='l2',
        datatype='',
        varformat=None,
        prefix='',
        suffix='',
        get_support_data=False,
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
    """
    This function loads SCM data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for SCM include ['brst' 'fast' 'slow' 'srvy']. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for SCM are: ['scsrvy', 'cal', 'scb', 'scf', 'schb', 'scm', 'scs']
            If no value is given the default is scsrvy for srvy data, and scb for brst data.

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)

            
    Returns:
        List of tplot variables created.

    """
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='scm',
            descriptor=datatype, varformat=varformat, prefix=prefix, suffix=suffix, get_support_data=get_support_data,
            time_clip=time_clip, no_update=no_update)
    return data,metadata


#TODO: Finish rewrite of mms_load_feeps
def mms_load_feeps(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy', 
        level='l2',
        datatype='electron',
        varformat=None,
        get_support_data=True,
        prefix='',
        suffix='',
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False,
        no_flatfield_corrections=False,
        data_units=['count_rate', 'intensity'] ):
    """
    This function loads FEEPS data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for FEEPS include ['brst', 'srvy']. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for FEEPS are: 
                       L2, L1b: ['electron', 'ion']
                       L1a: ['electron-bottom', 'electron-top', 'ion-bottom', 'ion-top']

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to True.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)

    Returns:
        List of tplot variables created.

    """
    #was: tvars = ...
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='feeps',
            descriptor=datatype, varformat=varformat, get_support_data=get_support_data, prefix=prefix, suffix=suffix,
            time_clip=time_clip, no_update=no_update)

    
    # Returning raw data until FEEPS-related sub-functions get migrated.
    return data,metadata

    
#    probes = probe if isinstance(probe, list) else [probe]
#    data_rates = data_rate if isinstance(data_rate, list) else [data_rate]
#    levels = level if isinstance(level, list) else [level]
#    datatypes = datatype if isinstance(datatype, list) else [datatype]
#
#    probes = [str(p) for p in probes]
#
#    #TODO: Migrate FEEPS-related sub-functions
#    mms_feeps_correct_energies(probes, data_rate, level=level, suffix=suffix)
#
#    if not no_flatfield_corrections:
#        mms_feeps_flat_field_corrections(probes=probes, data_rate=data_rate, suffix=suffix)
#
#    for probe in probes:
#        for datatype in datatypes:
#           mms_feeps_remove_bad_data(probe=probe, data_rate=data_rate, datatype =datatype, level=level, suffix=suffix)
#           
#           for data_unit in data_units:
#               eyes = mms_feeps_active_eyes(trange, probe, data_rate, datatype, level)
#
#               mms_feeps_split_integral_ch(data_unit, datatype, probe, suffix=suffix, data_rate=data_rate, level=level, sensor_eyes=eyes)
#
#               mms_feeps_remove_sun(eyes, trange, probe=probe, descriptor=datatype, data_units=data_unit, data_rate=data_rate, level=level, suffix=suffix)
#
#               mms_feeps_omni(eyes, probe=probe, descriptor=datatype, data_units=data_unit, data_rate=data_rate, level=level, suffix=suffix)
#
#               mms_feeps_spin_avg(probe=probe, data_units=data_unit, descriptor=datatype, data_rate=data_rate, level=level, suffix=suffix)
#
#    return tvars


# Function migrated from pyspedas: pyspedas/mms/eis/mms_eis_omni.py
def mms_eis_omni(
        probe,
        data_dictionary,
        species='proton',
        datatype='extof',
        suffix='',
        data_units='flux',
        data_rate='srvy', ):
    """
    This function will calculate the omni-directional EIS spectrograms, and is automatically called from mms_load_eis
    
    Parameters:
        probe: str
            probe #, e.g., '4' for MMS4

        data_dictionary: dict
            Dictionary of variables to operate on.

        species: str
            species for calculation (default: 'proton')

        datatype: str
            'extof' or 'phxtof' (default: 'extof')

        suffix: str
            suffix of the loaded data

        data_units: str
            'flux' or 'cps' (default: 'flux')

        data_rate: str
            instrument data rate, e.g., 'srvy' or 'brst' (default: 'srvy')


    Returns:
        Name of tplot variable created.
    """
    
    probe = str(probe)
    species_str = datatype + '_' + species
    if data_rate == 'brst':
        prefix = 'mms' + probe + '_epd_eis_brst_'
    else: 
        prefix = 'mms' + probe + '_epd_eis_'

    if data_units == 'flux':
        units_label = '1/(cm^2-sr-s-keV)'
    elif data_units == 'cps':
        units_label = '1/s'
    elif data_units == 'counts':
        units_label = 'counts'

    # telescopes = tnames(pattern=prefix + species_str + '_*' + data_units + '_t?'+suffix)
    pattern = prefix + species_str + '_*' + data_units + '_t?'+suffix
    telescopes = list()
    telescopes.extend(fnmatch.filter(data_dictionary.keys, pattern))

    if len(telescopes) == 6:
        time, data, energies = data_dictionary[telescopes[0]].values()
        flux_omni = np.zeros((len(time), len(energies)))
        for t in telescopes:
            time, data, energies = data_dictionary[t].values()
            flux_omni += data

        data_dictionary[prefix + species_str + '_' + data_units + '_omni' + suffix] = {'x': time, 'y': flux_omni/6., 'v': energies}
        # Removing redundant return value.
        #return prefix + species_str + '_' + data_units + '_omni' + suffix
    else:
        print('Error, problem finding the telescopes to calculate omni-directional spectrograms')
        # Removing redundant return value.
        #return None


def mms_load_eis(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy',
        level='l2',
        datatype='phxtof',
        varformat=None,
        get_support_data=False,
        prefix='',
        suffix='',
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
    """
    This function loads EIS data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for EIS include ['brst', 'srvy']. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for EIS are: ['extof', 'phxtof', and 'electronenergy']; default is 'extof'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)
            
    Returns:
        List of tplot variables created.

    """

    from .eis.mms_eis_omni import mms_eis_omni
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='epd-eis',
            descriptor=datatype, varformat=varformat, get_support_data=get_support_data, prefix='', suffix='',
            time_clip=time_clip, no_update=no_update)

    if data == []:
        return data,metadata

    if not isinstance(probe, list): probe = [probe]
    if not isinstance(data_rate, list): data_rate = [data_rate]
    if not isinstance(datatype, list): datatype = [datatype]

    for probe_id in probe:
        for datatype_id in datatype:
            for data_rate_id in data_rate:
                mms_eis_omni(probe_id, data_dictionary=data, data_rate=data_rate_id, datatype=datatype_id)
    return data,metadata



def mms_load_edi(
        trange=['2016-10-16', '2016-10-17'],
        probe='1',
        data_rate='srvy',
        level='l2',
        datatype='efield',
        varformat=None,
        get_support_data=False,
        prefix='',
        suffix='',
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
    """
    This function loads EDI data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for EDI include ['brst', 'fast', 'slow', 'srvy']. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for EDI are: ['efield', 'amb']; default is 'efield'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)
            
    Returns:
        List of tplot variables created.

    """
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='edi',
            descriptor=datatype, varformat=varformat, get_support_data=get_support_data, prefix=prefix, suffix=suffix, time_clip=time_clip, no_update=no_update)
    return data,metadata



def mms_load_edp(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='fast',
        level='l2',
        datatype='dce',
        varformat=None,
        get_support_data=False,
        prefix='',
        suffix='',
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
    """
    This function loads EDP data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for EDP include ['brst', 'fast', 'slow', 'srvy']. The
            default is 'fast'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for EDP are: ['dce', 'dcv', 'ace', 'hmfe']; default is 'dce'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)
            
    Returns:
        List of tplot variables created.

    """
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='edp',
            descriptor=datatype, varformat=varformat, get_support_data=get_support_data, prefix=prefix, suffix=suffix,
            time_clip=time_clip, no_update=no_update)
    return data,metadata



def mms_load_dsp(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy', 
        level='l2',
        datatype='bpsd',
        varformat=None,
        prefix='',
        suffix='',
        get_support_data=False,
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
    """
    This function loads DSP data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for DSP include ['fast', 'slow', 'srvy']. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for DSP are: ['epsd', 'bpsd', 'swd']; default is 'bpsd'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)
            
    Returns:
        List of tplot variables created.

    """
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='dsp',
            descriptor=datatype, varformat=varformat, prefix=prefix, suffix=suffix, get_support_data=get_support_data, time_clip=time_clip, no_update=no_update)
    return data,metadata



def mms_load_aspoc(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='srvy', 
        level='l2',
        datatype='',
        varformat=None,
        get_support_data=False,
        prefix='',
        suffix='',
        time_clip=False,
        no_update=False,
        available=False,
        notplot=False ):
    """
    This function loads ASPOC data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            instrument data rates for ASPOC include 'srvy', 'sitl'. The
            default is 'srvy'.

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatypes for ASPOC are: ['asp1', 'asp2', 'aspoc']; default is 'aspoc'

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

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

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)

    Returns:
        List of tplot variables created.

    """
    data,metadata = mms_load_data(trange=trange, notplot=notplot, probe=probe, data_rate=data_rate, level=level, instrument='aspoc',
            descriptor=datatype, varformat=varformat, get_support_data=get_support_data, prefix=prefix, suffix=suffix,
            time_clip=time_clip, no_update=no_update)
    return data,metadata



def mms_load_fsm(
        trange=['2015-10-16', '2015-10-17'],
        probe='1',
        data_rate='brst', 
        level='l3',
        datatype='8khz',
        get_support_data=False,
        time_clip=False,
        no_update=False, 
        available=False,
        varformat=None,
        notplot=False ):
    """
    This function loads FSM data into tplot variables
    
    Parameters:
        trange : list of str
            time range of interest [starttime, endtime] with the format 
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day 
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        probe : str or list of str
            list of probes, valid values for MMS probes are ['1','2','3','4']. 

        data_rate : str or list of str
            the current instrument data rate for FSM is 'brst'

        level : str
            indicates level of data processing. the default if no level is specified is 'l2'

        datatype : str or list of str
            Valid datatype for FSM is: 8khz

        get_support_data: bool
            If True, data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into data tables.  If False, only loads in data with a 
            "VAR_TYPE" attribute of "data".  Defaults to False.

        time_clip: bool
            Data will be clipped to the exact trange specified by the trange keyword.
            
        prefix: str
            The tplot variable names will be given this prefix.  By default, 
            no prefix is added.

        suffix: str
            The tplot variable names will be given this suffix.  By default, 
            no suffix is added.

        notplot: bool
            If True, then data are returned in a hash table instead of 
            being stored in tplot variables (useful for debugging, and
            access to multi-dimensional data products)
            
    Returns:
        List of tplot variables created.

    """
    data,metadata = mms_load_data(trange=trange, notplot=notplot, varformat=varformat, probe=probe, data_rate=data_rate, level=level, instrument='fsm', descriptor=datatype, get_support_data=get_support_data, time_clip=time_clip, no_update=no_update)
    return data,metadata
