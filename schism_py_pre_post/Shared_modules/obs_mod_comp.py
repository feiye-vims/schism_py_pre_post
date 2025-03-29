import pandas as pd
from scipy import interpolate
import numpy as np
from datetime import datetime, timezone
import copy


class obs_mod_comp():
    """
    class for handling model-data comparision
    """
    def __init__(self, obs=None, mod=None):
        self.obs_df = obs
        self.mod_df = mod
        self.obs_dt = None
        self.mod_dt = None
        self.mod_interp_df = None
        self.mod_interp_dt = None
        self.valid = False

        tmp_dfs = [self.obs_df, self.mod_df]
        tmps = [obs, mod]
        for i, _ in enumerate(tmps):
            tmp_dfs[i].rename(columns={tmp_dfs[i].columns[0]: "time"}, inplace=True)
            tmp_dfs[i].rename(columns={tmp_dfs[i].columns[1]: "value"}, inplace=True)

        obs_day = self.dattime2day(self.obs_df["time"], start_date=self.mod_df['time'][0])
        self.obs_dt = obs_day[1] - obs_day[0]
        self.obs_df.insert(1, 'Day', obs_day, True)

        mod_day = self.dattime2day(self.mod_df["time"], start_date=self.mod_df['time'][0])
        self.mod_dt = mod_day[1] - mod_day[0]
        self.mod_df.insert(1, 'Day', mod_day, True)

        obs_in_mod_time = (self.obs_df['Day'].to_numpy() >= self.mod_df['Day'].to_numpy()[0]) & \
                          (self.obs_df['Day'].to_numpy() <= self.mod_df['Day'].to_numpy()[-1])

        # sometimes no obs data is available in the model time range
        if obs_in_mod_time.sum() == 0:
            raise ValueError("No observation data available in the model time range")

        self.obs_df = self.obs_df[obs_in_mod_time]

        self.mod_interp_df = self.interp_mod_to_obs_time()
        mod_interp_day = self.mod_interp_df["time"].values
        self.mod_interp_dt = mod_interp_day[1] - mod_interp_day[0]
        self.mod_interp_df.insert(1, 'Day', mod_interp_day, True)

        self.valid = True

    def get_moving_average(self, nday_avg=3):
        self.obs_df['value'].values[:] = self.moving_average(self.obs_df['value'].values, n=int(nday_avg/self.obs_dt))
        self.mod_df['value'].values[:] = self.moving_average(self.mod_df['value'].values, n=int(nday_avg/self.mod_dt))
        self.mod_interp_df['value'].values[:] = self.moving_average(self.mod_interp_df['value'].values, n=int(nday_avg/self.mod_interp_dt))
        pass

    @staticmethod
    def dattime2day(dates, start_date=None, name='Day'):
        """
        Converts a series of dates to a series of float values that represent days since start_date.
        """
        if start_date:
            ts0 = pd.Timestamp(start_date).timestamp()
        else:
            ts0 = 0
        return ((dates.apply(pd.Timestamp.timestamp) - ts0)/(24*3600)).rename(name)

    @staticmethod
    def moving_average(a, n=3):
        ret = np.cumsum(a, axis=0, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ret[n-1:] = ret[n-1:] / n

        # re-align time series
        ret1 = ret * 0.0
        m = int(np.floor(n/2))
        ret1[m:-m] = ret[2*m:]

        # fill the first and last few records
        ret1[:m] = ret1[m]
        ret1[-m:] = ret1[-m-1]

        return ret1

    def interp_mod_to_obs_time(self):
        f_interp = interpolate.interp1d(self.mod_df['Day'], self.mod_df['value'])
        mod_interp_df = pd.DataFrame({'time': self.obs_df['Day'], 'value': f_interp(self.obs_df['Day'])})
        return mod_interp_df

    def cal_stats(self):
        import sklearn.metrics as metrics  # scikit-learn

        i_valid = np.logical_not((np.isnan(self.mod_interp_df['value'].to_numpy())) | np.isnan(self.obs_df['value'].to_numpy()))
        if np.any(i_valid):
            yhat = self.mod_interp_df['value'][i_valid].to_numpy()  # prediction
            y = self.obs_df['value'][i_valid].to_numpy()  # observation

            yhat_demeaned = yhat - np.mean(yhat)
            y_demeaned = y - np.mean(y)

            d = yhat - y

            mae = metrics.mean_absolute_error(y, yhat)

            mse = metrics.mean_squared_error(y, yhat)
            rmse = np.sqrt(mse)

            un_biased_mse = metrics.mean_squared_error(y_demeaned, yhat_demeaned)
            unbiased_rmse = np.sqrt(un_biased_mse)

            CC = np.corrcoef(y, yhat)[0, 1]

            stats_dict = {
                'Bias': np.mean(d),
                'MAE': mae,
                'RMSE': rmse,
                'CC': CC,
                'ubRMSE': unbiased_rmse,
            }
        else:
            stats_dict = {
                'Bias': np.nan,
                'MAE': np.nan,
                'RMSE': np.nan,
                'CC': np.nan,
                'ubRMSE': np.nan,
            }

        self.stats_dict = stats_dict
        self.stats_str = copy.copy(stats_dict)
        for key in stats_dict:
            self.stats_str[key] = "{:.3f}".format(stats_dict[key]),

        return stats_dict


if __name__ == "__main__":
    obs_times = [x.replace(tzinfo=timezone.utc) for x in pd.date_range('2021-03-21', periods=45, freq='1D')]
    obs_data = np.linspace(0, 2, 45)
    obs_df = pd.DataFrame({'datetime': obs_times, 'value': obs_data})

    mod_times = [x.replace(tzinfo=timezone.utc) for x in pd.date_range('2021-04-01', periods=15*24, freq='1H')]
    mod_data = np.ones(15*24)
    mod_df = pd.DataFrame({'datetime': mod_times, 'value': mod_data})

    my_comp = obs_mod_comp(obs=obs_df, mod=mod_df)
    stats_dict = my_comp.cal_stats()

    pass
