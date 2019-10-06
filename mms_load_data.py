#!/usr/bin/python
# -*- coding: utf-8 -*-

# This file adapted from mms_load_data.py from the pyspedas library,
# sourced from https://github.com/spedas/pyspedas
#
# All modifications copyright 2019 Tim Rogers.  All rights reserved.
# Released under the MIT license.
#
# Original copyright notice from pyspedas source preserved below:

# Copyright (c) 2017, THEMIS group, Space Sciences Laboratory, UC Berkeley.
# 
# All rights reserved.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import os
import requests
import logging
import numpy as np
from load_cdf import load_cdf
import mms_load_data_shims
from p_tqdm import p_map

#from ..spdtplot.cdf_to_tplot import cdf_to_tplot
#from ..analysis.time_clip import time_clip as tclip
#from pyspedas import time_double, time_string
from dateutil.parser import parse
from datetime import timedelta, datetime, timezone
#from shutil import copyfileobj, copy
#from tempfile import NamedTemporaryFile
#from .mms_config import CONFIG
#from .mms_get_local_files import mms_get_local_files
#from .mms_files_in_interval import mms_files_in_interval
#from .mms_login_lasp import mms_login_lasp

import mms_sdc_api_client

logging.captureWarnings(True)
logging.basicConfig(format='%(asctime)s: %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def mms_load_data(trange=['2015-10-16', '2015-10-17'], probe='1', data_rate='srvy', level='l2', 
    instrument='fgm', datatype='', descriptor=None, varformat=None, prefix='', suffix='', get_support_data=False, time_clip=False, 
    no_update=False, center_measurement=False, notplot=False, data_root=None):
    """
    This function loads MMS data into a dictionary by variable name.
    """

    if not isinstance(probe, list): probe = [probe]
    if not isinstance(data_rate, list): data_rate = [data_rate]
    if not isinstance(level, list): level = [level]
    if not isinstance(datatype, list): datatype = [datatype]
    if not isinstance(descriptor, list): descriptor = [descriptor]
    
    probe = [('mms'+(str(p))) for p in probe]

    ## allows the user to pass in trange as list of datetime objects
    #if type(trange[0]) == datetime and type(trange[1]) == datetime:
    #    trange = [time_string(trange[0].timestamp()), time_string(trange[1].timestamp())]
    #    
    #start_date = parse(trange[0]).strftime('%Y-%m-%d') # need to request full day, then parse out later
    #end_date = parse(time_string(time_double(trange[1])-0.1)).strftime('%Y-%m-%d-%H-%M-%S') # -1 second to avoid getting data for the next day
    
    # We're going to handle everything as datetime objects fo consistancy and easy conversion at-need.
    
    if type(trange[0]) == datetime: # Already a datetime.
        pass
    elif type(trange[0]) in (int,float,np.float32,np.float64): # Convert from posix timestamp if provided.
        trange[0] = datetime.fromtimestamp(trange[0], timezone.utc)
    elif type(trange[0]) == str: # Parse the string and generate a datetime.
        trange[0] = parse(trange[0])
    else:
        raise TypeError("Unsupported input format for start date/time.")
    
    if type(trange[1]) == datetime: # Already a datetime.
        pass
    elif type(trange[1]) == int: # Convert from posix timestamp if provided.
        trange[1] = datetime.fromtimestamp(trange[1], timezone.utc)
    elif type(trange[1]) == str: # Parse the string and generate a datetime.
        trange[1] = parse(trange[1])
    else:
        raise TypeError("Unsupported input format for end date/time.")
    
    # Replicating behavior of pyspedas:
    start_date = trange[0].date().isoformat() # need to request full day, then parse out later
    end_date = (trange[1] - timedelta(seconds=1)).isoformat() # -1 second to avoid getting data for the next day
    
    #download_only = CONFIG['download_only']
    #
    #no_download = False
    #if no_update or CONFIG['no_download']: no_download = True
    #
    #user = None
    #if not no_download:
    #    sdc_session, user = mms_login_lasp()
    
    out_files = []
    #available_files = []
    
    for dtype in datatype:
        # Default to 'science' data, as the old SPEDAS implementation assumed that and used "datatype" for something else entirely.
        if len(dtype) == 0: dtype = 'science'
        for lvl in level:
            for desc in descriptor:
                mms_api_client = mms_sdc_api_client.MMS_SDC_API_CLIENT(
                        sc=None, 
                        instr=instrument, 
                        mode=data_rate, 
                        level=lvl,
                        data_type=dtype,
                        data_root=data_root,
                        end_date=end_date,
                        offline=no_update,
                        optdesc=desc,
                        site='public',
                        start_date=start_date)
                logging.info('download URI: '+mms_api_client.url())
                out_files.extend(mms_api_client.Download())

    #for prb in probe:
    #    for drate in data_rate:
    #        for lvl in level:
    #            for dtype in datatype:
    #                if user is None:
    #                    url = 'https://lasp.colorado.edu/mms/sdc/public/files/api/v1/file_info/science?start_date=' + start_date + '&end_date=' + end_date + '&sc_id=mms' + prb + '&instrument_id=' + instrument + '&data_rate_mode=' + drate + '&data_level=' + lvl
    #                else:
    #                    url = 'https://lasp.colorado.edu/mms/sdc/sitl/files/api/v1/file_info/science?start_date=' + start_date + '&end_date=' + end_date + '&sc_id=mms' + prb + '&instrument_id=' + instrument + '&data_rate_mode=' + drate + '&data_level=' + lvl
    #                
    #                if dtype != '':
    #                    url = url + '&descriptor=' + dtype
    #
    #                if CONFIG['debug_mode']: logging.info('Fetching: ' + url)
    #
    #                if no_download == False:
    #                    # query list of available files
    #                    try:
    #                        http_json = sdc_session.get(url, verify=True).json()
    #
    #                        if CONFIG['debug_mode']: logging.info('Filtering the results down to your trange')
    #
    #                        files_in_interval = mms_files_in_interval(http_json['files'], trange)
    #
    #                        for file in files_in_interval:
    #                            file_date = parse(file['timetag'])
    #                            if dtype == '':
    #                                out_dir = os.sep.join([CONFIG['local_data_dir'], 'mms'+prb, instrument, drate, lvl, file_date.strftime('%Y'), file_date.strftime('%m')])
    #                            else:
    #                                out_dir = os.sep.join([CONFIG['local_data_dir'], 'mms'+prb, instrument, drate, lvl, dtype, file_date.strftime('%Y'), file_date.strftime('%m')])
    #
    #                            if drate.lower() == 'brst':
    #                                out_dir = os.sep.join([out_dir, file_date.strftime('%d')])
    #
    #                            out_file = os.sep.join([out_dir, file['file_name']])
    #
    #                            if CONFIG['debug_mode']: logging.info('File: ' + file['file_name'] + ' / ' + file['timetag'])
    #
    #                            if os.path.exists(out_file) and str(os.stat(out_file).st_size) == str(file['file_size']):
    #                                if not download_only: logging.info('Loading ' + out_file)
    #                                out_files.append(out_file)
    #                                continue
    #
    #                            if user is None:
    #                                download_url = 'https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?file=' + file['file_name']
    #                            else:
    #                                download_url = 'https://lasp.colorado.edu/mms/sdc/sitl/files/api/v1/download/science?file=' + file['file_name']
    #
    #                            logging.info('Downloading ' + file['file_name'] + ' to ' + out_dir)
    #
    #                            fsrc = sdc_session.get(download_url, stream=True, verify=True)
    #                            ftmp = NamedTemporaryFile(delete=False)
    #
    #                            with open(ftmp.name, 'wb') as f:
    #                                copyfileobj(fsrc.raw, f)
    #
    #                            if not os.path.exists(out_dir):
    #                                os.makedirs(out_dir)
    #
    #                            # if the download was successful, copy to data directory
    #                            copy(ftmp.name, out_file)
    #                            out_files.append(out_file)
    #                            fsrc.close()
    #                            ftmp.close()
    #                    except requests.exceptions.ConnectionError:
    #                        # No/bad internet connection; try loading the files locally
    #                        logging.error('No internet connection!')
    #
    #                
    #                if out_files == []:
    #                    if not download_only: logging.info('Searching for local files...')
    #                    out_files = mms_get_local_files(prb, instrument, drate, lvl, dtype, trange)
    #
    #if not no_download:
    #    sdc_session.close()
    #
    #if not download_only:
    out_files = sorted(out_files)
    
    #Because we're not using pytplot, new_variables is a simple dictionary containing the data.
    #  eg.
    #    pytplot:  get_data(Varname)
    #    current:  new_variables[Varname].values()
    #new_variables,new_metadata = load_cdf(out_files, varformat=varformat, get_support_data=get_support_data, prefix=prefix, suffix=suffix, center_measurement=center_measurement)
    
    ##alt:
    new_variables = {}
    new_metadata = {}
    
    logging.info('Beginning parallel load of '+str(len(out_files))+' data files...')
    
    pile_o_data = p_map(load_cdf, out_files, varformat, get_support_data, prefix, suffix, center_measurement)

    ########### Do merge
    logging.info('Stitching together the data...')
    for data,metadata in pile_o_data:
        # merge data dictionary
        for dataset in data.keys():
            if dataset not in new_variables.keys():
                # No previous entries for this dataset.  Add the whole thing as-is.
                new_variables[dataset] = data[dataset]
            else:
                # Panic and error out if the datasets with identical names don't have the same axes/variables being tracked.
                if len(new_variables[dataset].keys()) != len(data[dataset].keys()):
                    logging.error('Failure to merge new data with axes ('+(','.join((data[dataset].keys())))+') with existing data with axes ('+(','.join((new_variables[dataset].keys())))+')'+'.')
                    raise TypeError('Failure to merge new data with axes ('+(','.join((data[dataset].keys())))+') with existing data with axes ('+(','.join((new_variables[dataset].keys())))+')'+'.')
                
                # Update existing dataset entry with the additional data.
                for axis in data[dataset].keys():
                    new_variables[dataset][axis] = np.concatenate((new_variables[dataset][axis],data[dataset][axis]))
        # write/revise metadata
        for dataset in metadata.keys():
            if dataset not in new_metadata.keys():
                # No previous entries for this dataset.  Add the whole thing as-is.
                new_metadata[dataset] = metadata[dataset]
            else:
                # Compare the new set's metadata with the existing metadata.
                for meta in [key for key in metadata[dataset].keys() if key in set(new_metadata[dataset].keys())]:
                    # Print a notice for any unexpected differences, but don't raise exceptions.
                    if metadata[dataset][meta] != new_metadata[dataset][meta]:
                        logging.warning("Dataset '"+dataset+"' has non-matching metadata between input files. Old: {'"+meta+"': '"+new_metadata[dataset][meta]+"'} -- Now using: {'"+meta+"': '"+metadata[dataset][meta]+"'}")
                # Update the metadata, overwriting old values if appliciable.
                new_metadata[dataset].update(metadata[dataset])
    
    if len(new_variables) == 0:
        logging.warning('No data loaded.')
        return
    
    if len(new_metadata) == 0:
        logging.warning('No metadata loaded.')
        return
    
    logging.info('Loaded variables:')
    for new_var in new_variables.keys():
        print(new_var)
    
    if time_clip:
        logging.info('Clipping variables to time range...')
        mms_data_time_clip(new_variables, trange[0], trange[1])
    
    return new_variables, new_metadata
    #else:
    #    return out_files


def mms_data_time_clip(data_dict, start_time, end_time):
    """
    Revised function to clip all datasets in the data_dict to the specified time range.
    
    This function assumes the following:
        - in each dataset, there exists a numpy array 'x', which is a time axis
        - all values in the time axis are of the same type
        - being numpy arrays, the aforementioned type will be some form of float, representing a unix timestamp
    
    """

    new_start = 0.
    new_end = 0.

    # Handle likely input types for start and end times.  Convert to unix timestamp.
    if type(start_time) == datetime: # Handle all datetimes as UTC timezone
        new_start = start_time.replace(tzinfo=timezone.utc).timestamp()
    elif type(start_time)  in (int,float,np.float32,np.float64): # Looks like it's already a timestamp. Use as-is.
        new_start = start_time
    elif type(start_time) == str: # Parse the string, then convert to timestamp.
        new_start = parse(start_time).replace(tzinfo=timezone.utc).timestamp()
    else:
        raise TypeError("Unsupported input format for start date/time.")
    
    if type(end_time) == datetime: # Handle all datetimes as UTC timezone
        new_end = end_time.replace(tzinfo=timezone.utc).timestamp()
    elif type(end_time)  in (int,float,np.float32,np.float64): # Looks like it's already a timestamp. Use as-is.
        new_end = end_time
    elif type(end_time) == str: # Parse the string, then convert to timestamp.
        new_end = parse(end_time).replace(tzinfo=timezone.utc).timestamp()
    else:
        raise TypeError("Unsupported input format for end date/time.")
    
    
    for setname in data_dict.keys():
        dataset = data_dict[setname]
        timearr = dataset['x']
        start_index = 0
        end_index = len(timearr)-1
        
        # Early sanity checks
        if new_start > new_end:
            logging.error('Error: Start time is after end time.')
            continue
        if (new_start > timearr[-1]) or (new_end < timearr[0]):
            logging.warning('Entire dataset removed by time clipping: "'+setname+'"')
            continue
        if (new_start <= timearr[0]) and (new_end >= timearr[-1]):
            logging.info('Entire dataset included in time clip range: "'+setname+'"')
            continue
        
        # Find our indices
        while timearr[start_index] < new_start: start_index += 1
        while timearr[end_index] > new_end: end_index -= 1
        
        # End index needs to point at the first index we want to exclude after the data.
        end_index += 1

        # Update the datasets
        for axis in dataset.keys():
            # Only clip the axis if it is a list/array instead of a scalar.
            if len(dataset[axis].shape) >= 1:
                dataset[axis] = dataset[axis][start_index:end_index]

    # Return the updated data dictionary
    return data_dict
