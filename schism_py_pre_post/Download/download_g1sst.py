import os
import pandas as pd


def download_g1sst(
    download_dir='./', begin_date='2015-09-26', end_date='2015-09-28',
    box=[-99, -60, 7, 47],
    valid_size=109436146, timeout=500,
    i_extract=True
):
    '''
    sample url:
    https://www.ncei.noaa.gov/data/oceans/ghrsst/L4/GLOB/JPL_OUROCEAN/G1SST/2015/001/20150101-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2
    '''
    # convert to datetime
    days = pd.date_range(begin_date, end_date)

    for i in range(0, 3):  # multiple tries
        for this_day in days:
            outname = f'{download_dir}/{this_day.strftime("%Y%m%d")}.nc.bz2'
            print('---------------------------------------------------------------------')
            if not os.path.exists(outname) or os.stat(outname).st_size < valid_size*0.92:
                url = (
                       'https://www.ncei.noaa.gov/data/oceans/ghrsst/L4/GLOB/JPL_OUROCEAN/G1SST/'
                       + this_day.strftime("%Y") + '/' + str(this_day.day_of_year).zfill(3) + '/' + this_day.strftime("%Y%m%d")
                       + '-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2'
                      )
                print(url)
                os.system(f'curl -m {timeout}  "{url}" > {outname}')
                if i_extract and not os.path.exists(outname[:-4]):
                    print(f'extracting from {outname}')
                    os.system(f'bzip2 -dk {outname}')
                    
            else:
                print('file exists: ' + outname)

