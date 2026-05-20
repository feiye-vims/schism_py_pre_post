import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pylib_experimental.schism_file import TimeHistory
from schism_py_pre_post.Grid.Bpfile import Bpfile


def mad(x: np.ndarray) -> float:
    """Median Absolute Deviation (MAD)."""
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return np.median(np.abs(x - med))


def artificial_data(h: pd.Series) -> pd.Series:
    """Generate artificial data with spikes and plateaus for testing."""
    h_art = h.copy()

    h_art = h.loc["2022-09-25":"2022-10-01"].copy()

    plateau_start = pd.Timestamp("2022-09-26 06:00:00")
    plateau_end = pd.Timestamp("2022-09-26 12:00:00")

    h_art.loc[plateau_start: plateau_end] = -2  # plateau glitch

    return h_art


def mark_jumping_plateaus(
    h: pd.Series, smooth_window: int, max_gap_samples: int,
    jump_tol: float, plateau_side_diff_tol: float,
):
    """
    Identify jump-plateau-jump glitches in a time series.
    
    inputs:
        h : pd.Series, time series with DatetimeIndex
        smooth_window : int, rolling median window (in samples) to calculate per-sample change
        max_gap_samples : int, maximum number of samples between up-jump and down-jump to treat as plateau
        jump_tol : float, threshold on robust derivative to identify jumps
        plateau_side_diff_tol : float, maximum allowed difference between medians before/after plateau
    """

    h_med = h.rolling(window=smooth_window, center=True, min_periods=1).median()
    dh = h_med.diff()  # per-sample change

    up_jumps = dh > jump_tol
    down_jumps = dh < -jump_tol

    up_idx = np.where(up_jumps.values)[0]
    down_idx = np.where(down_jumps.values)[0]

    max_gap_samples = max(1, max_gap_samples)

    plateau_mask = np.zeros(len(h), dtype=bool)

    def _mark_plateau(i: int, j: int) -> None:
        """Mark [i, j] as plateau if it returns close to baseline."""
        if (j - i) > max_gap_samples:
            return
        before = h.iloc[max(i - 5, 0): i].median()
        after = h.iloc[j: min(j + 5, len(h))].median()
        if np.isfinite(before) and np.isfinite(after) and (np.abs(after - before) < plateau_side_diff_tol):
            # similar medians before/after jumps, i.e., plateau
            plateau_mask[i: j + 1] = True

    # (1) up -> down
    for i in up_idx:
        downs_after = down_idx[down_idx > i]
        if downs_after.size == 0:
            continue
        j = int(downs_after[0])
        _mark_plateau(i, j)

    # (2) down -> up
    for i in down_idx:
        ups_after = up_idx[up_idx > i]
        if ups_after.size == 0:
            continue
        j = int(ups_after[0])
        _mark_plateau(i, j)

    return plateau_mask


def mark_impulse_spikes(
    h: pd.Series, window: int, k_impulse: float,
) -> pd.Series:
    """
    Identify impulse spikes in a time series using Hampel-style detection.

    inputs:
        h : pd.Series, time series with DatetimeIndex
        window : int, rolling window (in samples) for local median/MAD
        k_impulse : float, multiplier on local robust scale for spike detection
    """

    rolling_med = h.rolling(window=window, center=True, min_periods=1).median()
    rolling_mad = h.rolling(window=window, center=True, min_periods=1).apply(lambda x: mad(x.values), raw=False)

    local_sigma = 1.4826 * rolling_mad
    # Avoid zeros
    local_sigma = local_sigma.where(local_sigma > 0, other=np.nan)

    impulse_mask = (np.abs(h - rolling_med) > (k_impulse * local_sigma)).fillna(False).values

    return impulse_mask


def remove_spikes_and_plateaus(
    h: pd.Series,
    smooth_window: int = 5, hampel_window: int = 7,
    k_jump: float = 6.0, k_impulse: float = 6.0,
    max_gap: str = "8h", return_tol: float = 2.0,
) -> tuple[pd.Series, pd.Series]:
    """
    Remove (set to NaN) both:
      1) impulse spikes (Hampel-style median/MAD outliers), and
      2) jump-plateau-jump glitches (paired opposite-signed steps).

    Parameters
    ----------
    h : pd.Series
        Time series with a DatetimeIndex.
    smooth_window : int
        Rolling median window (in samples) for calculating per-sample change.
    hampel_window : int
        Rolling window (in samples) for impulse spike detection.
    k_jump : float
        Multiplier on robust scale (sigma) for jump edge detection.
    k_impulse : float
        Multiplier on local robust scale for impulse spike detection.
    max_gap : str
        Maximum duration between an up-jump and the subsequent down-jump
        to treat as a plateau glitch (e.g., '2h', '30min').
    return_tol : float
        Require that the median after the down-jump returns close to the
        median before the up-jump: |after-before| < return_tol * sigma.

    Returns
    -------
    h_clean : pd.Series
        Series with detected bad points/intervals set to NaN.
    mask : pd.Series[bool]
        Boolean mask (True = removed).
    """
    if not isinstance(h.index, pd.DatetimeIndex):
        raise TypeError("h must have a DatetimeIndex")

    # Ensure monotonic time index
    h = h.sort_index()
    # Sampling interval (assumes mostly regular)
    if len(h.index) < 2:
        return h.copy(), pd.Series(False, index=h.index)

    dt = h.index[1] - h.index[0]
    if dt <= pd.Timedelta(0):
        raise ValueError("Non-positive time step detected")

    # Global robust scale (sigma-like)
    sigma = 1.4826 * mad(h.values)
    if not np.isfinite(sigma) or sigma == 0:
        # fallback: small epsilon to avoid division-by-zero
        sigma = np.nanstd(h.values)
        if not np.isfinite(sigma) or sigma == 0:
            sigma = 1e-12

    # --- (A) Jump / plateau detection via robust derivative ---
    plateau_mask = mark_jumping_plateaus(
        h,
        smooth_window=smooth_window,
        max_gap_samples=int(pd.Timedelta(max_gap) / dt),
        jump_tol=k_jump * sigma,
        plateau_side_diff_tol=return_tol * sigma,
    )

    # --- (B) Impulse spike detection (Hampel-style), currently not used ---
    impulse_mask = mark_impulse_spikes(
        h, window=hampel_window, k_impulse=k_impulse,
    )

    # --- (C) Jumps detection via derivative thresholding ---
    h_scaled = h / sigma
    dh_dt = h_scaled.diff() / (dt.total_seconds() / 3600.0)  # per-hour change
    dh2_dt2 = dh_dt.diff() / (dt.total_seconds() / 3600.0)  # 2nd derivative
    # spike_mask = (np.abs(dh_dt) > 0.5).fillna(False).values  # 0.5 m/h based on max observed jump rates
    spike_mask = (np.abs(dh2_dt2) > 1.8).fillna(False).values  # nan values (if any) will be treated as non-spikes

    # Combine
    mask = plateau_mask | spike_mask
    mask = pd.Series(mask, index=h.index)

    h_clean = h.copy()
    h_clean[mask] = np.nan

    # Fill nan with linear interpolation
    h_clean = h_clean.interpolate(method="time", limit_direction="both")

    return h_clean, mask


def plot_comparison(h, h_clean, mask, title_str: str):
    """
    Plot comparison of raw vs cleaned time series.
    """
    plt.figure()
    h.plot(label="raw", linewidth=1)
    h_clean.plot(label="clean", linewidth=1)

    # Optional: highlight removed points
    if mask.any():
        plt.scatter(h.index[mask], h[mask], s=8, marker="x", label="removed")

    plt.title(f"Station {title_str}: raw vs cleaned")
    plt.legend()
    plt.tight_layout()
    plt.show()


def single_station_test():
    station_in = "/sciclone/schism10/feiye/TEMP/EnOI/run01q/station.in"
    staout_1 = "/sciclone/schism10/feiye/STOFS3D-v7.3/r2022/staout_1"

    station_ids = Bpfile(station_in, cols=5).st_id

    th = TimeHistory.from_file(
        staout_1,
        start_time_str="2021-12-31",
        columns=station_ids,
    )

    # Pick station
    st = "8720226" # "8413320" 
    h = th.df[st]
    # h = artificial_data(h)

    # Clean
    h_clean, mask = remove_spikes_and_plateaus(
        h, smooth_window=5, hampel_window=7, k_jump=6.0, k_impulse=3.0,
        max_gap="12h", return_tol=2.0,
    )

    # Plot
    plot_comparison(h, h_clean, mask, title_str=st)


def main():
    station_in = "/sciclone/schism10/feiye/TEMP/EnOI/run01q/station.in"
    staout_1 = "/sciclone/schism10/feiye/STOFS3D-v7.3/Total/staout_1"
    
    station_ids = Bpfile(station_in, cols=5).st_id

    th = TimeHistory.from_file(staout_1, start_time_str="2021-12-31", columns=station_ids)

    for i, st in enumerate(station_ids):
        if i >= 164:
            break
        print(f"Processing station {st}...")
        h = th.df[st]
        h_clean, mask = remove_spikes_and_plateaus(
            h, smooth_window=5, hampel_window=7, k_jump=6.0, k_impulse=3.0,
            max_gap="12h", return_tol=2.0,
        )
        th.df[st] = h_clean
    th.writer(staout_1.replace("staout_1", "staout_1_cleaned"))


if __name__ == "__main__":
    # single_station_test()
    main()
    pass
