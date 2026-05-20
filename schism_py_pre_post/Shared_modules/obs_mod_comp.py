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
        self.obs_dt = None  # time step in the observation data, in days
        self.mod_dt = None  # time step in the model data, in days
        self.mod_interp_df = None
        self.mod_interp_dt = None
        self.valid = False

        tmp_dfs = [self.obs_df, self.mod_df]
        tmps = [obs, mod]
        for i, _ in enumerate(tmps):
            tmp_dfs[i].rename(columns={tmp_dfs[i].columns[0]: "time"}, inplace=True)
            tmp_dfs[i].rename(columns={tmp_dfs[i].columns[1]: "value"}, inplace=True)

        obs_seconds = self.obs_df["time"].apply(pd.Timestamp.timestamp) - self.mod_df["time"].iloc[0].timestamp()
        # round to whole seconds
        obs_seconds = np.round(obs_seconds).astype(int)
        obs_day = obs_seconds / (24*3600)
        self.obs_dt = self.infer_nominal_dt(obs_day)
        self.obs_df.insert(1, 'Day', obs_day, True)

        mod_seconds = self.mod_df["time"].apply(pd.Timestamp.timestamp) - self.mod_df["time"].iloc[0].timestamp()
        mod_day = mod_seconds / (24*3600)
        self.mod_dt = self.infer_nominal_dt(mod_day)
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
    
    def infer_nominal_dt(self, t_days):
        diffs = np.diff(t_days)
        dt0 = np.median(diffs)
        valid = diffs < 1.05 * dt0
        inferred_dt = np.mean(diffs[valid])
        return inferred_dt
    
    
    def make_clean_time_series(self, time_days, data, dt=None):
        """
        Map an irregular time series onto a regular grid with nominal time step.

        Behavior
        --------
        - Original observed samples are mapped to the nearest regular-grid index.
        - Original NaNs are preserved at their mapped positions.
        - Missing time steps remain NaN.
        - No interpolation is performed.

        Parameters
        ----------
        time_days : array-like
            Original time values in days.
        data : array-like
            Original data values.
        dt : float, optional
            Nominal time step in days. If None, infer from time_days.

        Returns
        -------
        time_clean : ndarray
            Regular time grid.
        data_clean : ndarray
            Data on the regular grid, with gaps and original NaNs retained.
        obs_mask : ndarray of bool
            True where an original sample was mapped onto the grid.
        dt : float
            Nominal time step used.
        """
        time_days = np.asarray(time_days, dtype=float)
        data = np.asarray(data, dtype=float)

        if len(time_days) != len(data):
            raise ValueError("time_days and data must have the same length.")
        if len(time_days) == 0:
            raise ValueError("time_days is empty.")

        if dt is None:
            dt = self.infer_nominal_dt(time_days)

        # Anchor grid at the first sample
        t0 = time_days[0]
        k = np.round((time_days - t0) / dt).astype(int)

        # Build full regular grid
        kmin, kmax = k.min(), k.max()
        time_clean = t0 + np.arange(kmin, kmax + 1) * dt
        data_clean = np.full(time_clean.shape, np.nan, dtype=float)
        obs_mask = np.zeros(time_clean.shape, dtype=bool)  # True where an original sample maps to the grid

        # Shift indices so they align with the local grid
        kk = k - kmin

        # Handle possible collisions: if multiple raw times map to same grid point,
        # keep the last finite value; if all are NaN, leave NaN.
        for i, j in enumerate(kk):
            obs_mask[j] = True
            if np.isfinite(data[i]):
                data_clean[j] = data[i]
            elif not np.isfinite(data_clean[j]):
                # preserve NaN if nothing finite has been assigned there
                data_clean[j] = np.nan

        return time_clean, data_clean, obs_mask, dt
    
    def split_into_segments(self, values):
        """
        Split indices into continuous segments based on NaNs in the data.

        Parameters
        ----------
        values : array-like
            Data values (NaNs indicate gaps).

        Returns
        -------
        segments : list of ndarray
            List of index arrays, one per continuous segment (no NaNs inside).
        """
        values = np.asarray(values, dtype=float)
        n = len(values)

        valid = ~np.isnan(values)

        segments = []
        i = 0

        while i < n:
            if not valid[i]:
                i += 1
                continue

            start = i
            while i < n and valid[i]:
                i += 1
            end = i

            segments.append(np.arange(start, end))

        return segments

    def lowpass_filter_segments(
        self,
        time_days,
        data,
        cutoff_period_days=2.0,
        filter_order=4,
        gap_factor=1.5,
        min_len_factor=3,
    ):
        """
        Apply zero-phase Butterworth low-pass filter to each continuous segment.

        Parameters
        ----------
        time_days : array-like
            Time values in days.
        data : array-like
            Data values.
        cutoff_period_days : float, optional
            Cutoff period in days. Default is 2 days.
        filter_order : int, optional
            Butterworth filter order.
        gap_factor : float, optional
            Gaps are defined as time jumps > gap_factor * dt.
        min_len_factor : int, optional
            Segment must have at least min_len_factor * padlen points.

        NaNs are linearly interpolated within the segment before filtering,
        and the filtered result recovers the NaN values at their original positions.
        

        Returns
        -------
        filtered : ndarray
            Filtered data, with NaNs in skipped regions.
        dt : float
            Inferred nominal time step in days.
        segments : list of ndarray
            Continuous segments used.
        """

        from scipy.signal import butter, filtfilt
        
        time_days = np.asarray(time_days, dtype=float)
        data = np.asarray(data, dtype=float)

        if len(time_days) != len(data):
            raise ValueError("time_days and data must have the same length.")

        time_clean, data_clean, obs_mask, dt = self.make_clean_time_series(time_days, data)
        segments = self.split_into_segments(data_clean)

        fs = 1.0 / dt  # samples per day
        fc = 1.0 / cutoff_period_days  # cycles per day

        b, a = butter(filter_order, fc, btype="low", fs=fs)

        # filtfilt default padlen is 3 * max(len(a), len(b))
        padlen = 3 * max(len(a), len(b))
        min_seg_len = max(padlen + 1, min_len_factor * padlen)

        filtered = np.full_like(data_clean, np.nan, dtype=float)

        for seg in segments:
            if len(seg) < min_seg_len:
                continue

            y = data_clean[seg].astype(float)

            if np.any(~np.isfinite(y)):
                raise ValueError("Unexpected: segment identified as valid but contains NaNs")
            
            if len(y) < min_seg_len:
                continue

            try:
                filtered[seg] = filtfilt(b, a, y)
            except ValueError:
                # Segment still too short or otherwise unsuitable
                continue

        return time_clean, filtered

    def apply_low_pass_filter(
        self,
        cutoff_period_days=2.0,
        filter_order=4,
        gap_factor=1.5,
        preserve_nan=False,
    ):
        """
        Apply low-pass filter with segment handling to obs/mod/mod_interp dataframes.

        Assumes:
        - self.obs_df, self.mod_df, self.mod_interp_df each contain:
            - a time column named 'day'
            - a value column named 'value'
        - time is in days

        Stores inferred dts into:
        - self.obs_dt
        - self.mod_dt
        - self.mod_interp_dt
        """
        # Observation
        time_clean, filtered = self.lowpass_filter_segments(
            time_days=self.obs_df["Day"].values,
            data=self.obs_df["value"].values,
            cutoff_period_days=cutoff_period_days,
            filter_order=filter_order,
            gap_factor=gap_factor,
        )
        self.obs_df["value"] = np.interp(self.obs_df["Day"].values, time_clean, filtered)

        # Model
        time_clean, filtered = self.lowpass_filter_segments(
            time_days=self.mod_df["Day"].values,
            data=self.mod_df["value"].values,
            cutoff_period_days=cutoff_period_days,
            filter_order=filter_order,
            gap_factor=gap_factor,
        )
        self.mod_df["value"] = np.interp(self.mod_df["Day"].values, time_clean, filtered)

        # Interpolated model
        time_clean, filtered = self.lowpass_filter_segments(
            time_days=self.mod_interp_df["Day"].values,
            data=self.mod_interp_df["value"].values,
            cutoff_period_days=cutoff_period_days,
            filter_order=filter_order,
            gap_factor=gap_factor,
        )
        self.mod_interp_df["value"] = np.interp(self.mod_interp_df["Day"].values, time_clean, filtered)

    def get_moving_average(self, nday_avg=3):
        self.obs_df['value'].values[:] = self.moving_average(self.obs_df['value'].values, n=int(nday_avg/self.obs_dt))
        self.mod_df['value'].values[:] = self.moving_average(self.mod_df['value'].values, n=int(nday_avg/self.mod_dt))
        self.mod_interp_df['value'].values[:] = self.moving_average(self.mod_interp_df['value'].values, n=int(nday_avg/self.mod_interp_dt))

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

            if np.std(y) == 0 or np.std(yhat) == 0:
                CC = np.nan
            else:
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
            self.stats_str[key] = "{:.3f}".format(stats_dict[key])

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
