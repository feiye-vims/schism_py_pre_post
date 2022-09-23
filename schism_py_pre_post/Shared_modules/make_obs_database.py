from schism_pre_post_process import ObsData
import time

if __name__ == "__main__":
    states = ['ME', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'DE', 'PA', 'MD', 'VA', 'NC', 'SC', 'GA', 'FL', 'AL', 'MI', 'LA', 'TX']
    states = ['AL', 'MI', 'LA', 'TX']
    for state in states:
        fname = f'/sciclone/schism10/feiye/schism20/REPO/USGS/25_yr_iv_discharge_1993_2018/{state}/{state}_discharge.pkl'

        print(state)
        t1 = time.time()
        my_obs = ObsData(fname, from_raw_data=False)
        print(time.time()-t1)
        pass

    pass
