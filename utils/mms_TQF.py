'''
Function(s) to extract Tetrahedron Quality Factor(TQF) and separation scale data
from local copies of MMS SDC files.  
'''

import pandas as pd
import numpy as np
import datetime as dt

def get_TQF(trange=['2017-05-01', '2017-10-01'], datadir="~/data/mms/ancillary/mms/defq/"):
    '''
    Extracts timestamp, TQF, and scale scieze for the tetrahedron (in km) from 
    locally stored files from the SDC.  This is meant to eventually pair with
    another function which will do appropriate version checking using the SDC
    API and download files as necessary.
    '''

    print('trange = ', trange)

    # First create a datetime object of the beginning time and stop time from trange
    btime = dt.datetime.strptime(trange[0][:10], "%Y-%m-%d")
    stime = dt.datetime.strptime(trange[1][:10], "%Y-%m-%d")
    if btime == stime: stime += dt.timedelta(days=1)

    # Convert this into string with format used in TQF filenames of YearDay
    btime = btime.strftime("%Y%j")
    stime = stime.strftime("%Y%j")

    # initialize list to hold filenames
    flist = []

    # initialize first date for filename
    t1 = int(btime)

    # loop to generate sequential filenames for days from beginning to end of trange
    while int(t1) < int(stime):
        if ((int(str(t1)[-3:]) + 1)%365) == (int(str(t1)[-3:])+1): t2 = int(t1)+1  # generates following day, accounting for end of year
        else: t2 = (int(str(t1)[:4])+1)*1000 +1 # NOTE: Does not correct for leap year
        flist.append("MMS_DEFQ_"+str(t1)+"_"+str(t2)+".V00")    # Does not correct for higher versions (future add?)
        t1 = int(t2)
    # endwhile

    # initialize pandas dataframe to hold loaded TQF data
    TQFdf = pd.DataFrame(columns=['Epoch(UTC)', 'TQF', 'Scale(km)']) # NOTE: Epoch is a data column here, NOT the index! (future change?)

    # step through the file list and append the data from each to the dataframe
    for fname in flist:
        file = pd.read_csv(datadir+fname, delim_whitespace=True, skiprows=11, names=['Epoch(UTC)', 'TQF', 'Scale(km)'], usecols=[0,2,3])
        file['Epoch(UTC)'] = file['Epoch(UTC)'].map(lambda date: dt.datetime.strptime(date[:-4], '%Y-%j/%H:%M:%S').timestamp() + float(date[-4:]))
        TQFdf = TQFdf.append(file)
    # endfor
    
    # clear duplicate entries from conceutive files
    TQFdf = TQFdf.drop_duplicates()

    # return both the dataframe and a numpy array; which is needed can be selected when function is called
    return (TQFdf, TQFdf.to_numpy())



