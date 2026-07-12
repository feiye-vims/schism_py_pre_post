"""
Plot 2D transect data, e.g., salinity along a river transect
"""
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from add_river_mile import RiverMilesRef
from pylib import read


def plot_transect_filled(
        x, z, var, var_name, t_idx=0,
        nx_i=600, nz_i=400,
        ax=None,
        cmap='viridis',
        reverse_x=False,
        draw_mesh=True,
        draw_color_bar=True,
        mesh_color='white',
        clim=None,
        mesh_alpha=0.6,
        mesh_lw=0.4,
        plot_args=None
    ):
    """
    Smoothly filled vertical transect at time index t_idx,
    masked above the surface AND below the bottom,
    with optional mesh lines.

    Works even if x is decreasing (e.g., river miles).
    """

    x = np.asarray(x)

    # --- extract time slice ---
    z2d   = z[t_idx]       # (nx, nz)
    z2d[z2d > 9999] = np.nan  # mask invalid depths
    var2d = var[t_idx]     # (nx, nz)
    nx, nz = z2d.shape

    # --- original coordinates (for plotting & griddata) ---
    X2d = np.repeat(x[:, None], nz, axis=1)  # (nx, nz)

    xs = X2d.ravel()
    zs = z2d.ravel()
    vs = var2d.ravel()

    # --- build uniform grid in physical x,z space ---
    x_i = np.linspace(xs.min(), xs.max(), nx_i)
    z_i = np.linspace(np.nanmin(zs), np.nanmax(zs), nz_i)
    Xi, Zi = np.meshgrid(x_i, z_i)

    # --- interpolate variable onto uniform grid ---
    Vi = griddata((xs, zs), vs, (Xi, Zi),
                  method='linear', fill_value=np.nan)

    # -----------------------
    #   CREATE MASKS
    # -----------------------

    # bottom(x) = deepest z at each x; surface(x) = shallowest
    bottom  = z2d.max(axis=1)   # (nx,)
    surface = z2d.min(axis=1)   # (nx,)

    # For interpolation we need x to be increasing.
    sort_idx = np.argsort(x)
    x_sorted       = x[sort_idx]
    bottom_sorted  = bottom[sort_idx]
    surface_sorted = surface[sort_idx]

    # interpolate bottom/surface onto the Xi grid in x
    bottom_i  = np.interp(x_i, x_sorted, bottom_sorted)     # (nx_i,)
    surface_i = np.interp(x_i, x_sorted, surface_sorted)    # (nx_i,)

    # mask deep region: z > bottom(x) and air region: z < surface(x)
    mask_below = Zi > bottom_i[np.newaxis, :]
    mask_above = Zi < surface_i[np.newaxis, :]
    mask = mask_below | mask_above

    Vi_masked = np.where(mask, np.nan, Vi)

    # -----------------------
    #     PLOTTING
    # -----------------------
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 5))

    if clim is None:
        vmin = np.nanmin(Vi_masked)
        vmax = np.nanmax(Vi_masked)
    else:
        vmin, vmax = clim
    pc = ax.pcolormesh(Xi, Zi, Vi_masked, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax)
    if draw_color_bar:
        cb = plt.colorbar(pc, ax=ax)
        cb.set_label(var_name)

    # mesh lines on top
    if draw_mesh:
        for j in range(nz):
            ax.plot(X2d[:, j], z2d[:, j],
                    color=mesh_color, alpha=mesh_alpha, linewidth=mesh_lw)
        for i in range(nx):
            ax.plot(X2d[i, :], z2d[i, :],
                    color=mesh_color, alpha=mesh_alpha, linewidth=mesh_lw)

    if reverse_x:
        ax.invert_xaxis()

    if plot_args is not None:
        ax.set(**plot_args)
    ax.set_xlabel("along-transect distance (miles)")
    ax.set_ylabel("bathymetry (m)")
    ax.set_title(f"Transect at t={t_idx}")
    ax.invert_yaxis()  # z positive downward

    # If user gave x in decreasing order (river miles), keep that orientation:
    if np.any(np.diff(x) < 0):
        ax.invert_xaxis()
    
    return ax


import numpy as np
from scipy.ndimage import uniform_filter1d

def cal_intrusion_length_from_mouth(
        x_native,
        salinity_native,
        threshold=0.5,
        dx_reg=250.0,
        smooth_km=1.0,
        min_consecutive=10,
        temporal_smooth_hours=12,
        dt_hours=1,
        return_front_x=False,
    ):
    """
    Calculate salt intrusion length along a transect.

    Definition
    ----------
    x_native must be ordered from upstream/fresh to downstream/mouth/salty:

        x_native[0]  = upstream end
        x_native[-1] = mouth / downstream end

    The returned intrusion length is measured from the mouth upstream:

        intrusion_length = x_mouth - x_front

    where x_front is the accepted upstream salinity front, defined by the
    threshold salinity.

    Therefore:
        - no salt anywhere       -> intrusion_length = 0
        - entire transect salty  -> intrusion_length = x_mouth - x_upstream

    Parameters
    ----------
    x_native : 1D array, shape (npoints,)
        Cumulative distance along transect in meters, increasing from
        upstream to downstream/mouth.

    salinity_native : 2D array, shape (nt, npoints)
        Bottom salinity along the transect.

    threshold : float
        Salinity threshold defining salt presence, e.g. 0.5 psu or 4 psu.

    dx_reg : float
        Regular-grid spacing in meters.

    smooth_km : float
        Spatial smoothing window in km.

    min_consecutive : int
        Minimum number of consecutive fresh points required upstream of the
        front. This suppresses isolated fresh/salty noise.

    temporal_smooth_hours : float or None
        If not None, apply running mean to intrusion_length.

    dt_hours : float or None
        Model output interval in hours.

    return_front_x : bool
        If True, also return the front location x_front.

    Returns
    -------
    intrusion_length : 1D array, shape (nt,)
        Salt intrusion length measured from the mouth, in same units as x_native.

    front_x : 1D array, shape (nt,), optional
        Accepted salinity-front location in x coordinates.
    """

    x_native = np.asarray(x_native)
    salinity_native = np.asarray(salinity_native)

    if salinity_native.ndim != 2:
        raise ValueError("salinity_native must be 2D with shape (nt, npoints).")

    nt, npoints = salinity_native.shape

    if x_native.ndim != 1:
        raise ValueError("x_native must be 1D.")

    if x_native.size != npoints:
        raise ValueError("x_native.size must match salinity_native.shape[1].")

    if not np.all(np.diff(x_native) > 0):
        raise ValueError(
            "x_native must be strictly increasing from upstream/fresh "
            "to downstream/mouth/salty."
        )

    # ---------------------------------------------------------
    # 1. Construct regular x-grid
    # ---------------------------------------------------------
    x_reg = np.arange(x_native[0], x_native[-1] + dx_reg, dx_reg)
    x_mouth = x_reg[-1]
    x_upstream = x_reg[0]
    n_reg = x_reg.size

    # ---------------------------------------------------------
    # 2. Interpolate salinity onto regular x-grid
    # ---------------------------------------------------------
    sal_reg = np.empty((nt, n_reg), dtype=float)

    for t in range(nt):
        sal_reg[t, :] = np.interp(x_reg, x_native, salinity_native[t, :])

    # ---------------------------------------------------------
    # 3. Spatial smoothing along x
    # ---------------------------------------------------------
    window_m = smooth_km * 1000.0
    window_pts = max(1, int(round(window_m / dx_reg)))

    sal_smooth = uniform_filter1d(
        sal_reg,
        size=window_pts,
        axis=1,
        mode="nearest",
    )

    # ---------------------------------------------------------
    # 4. Detect front and compute intrusion length from mouth
    # ---------------------------------------------------------
    front_x = np.full(nt, np.nan)
    intrusion_length = np.full(nt, np.nan)

    for t in range(nt):
        profile = sal_smooth[t, :]
        salty = profile >= threshold

        # No salt anywhere: intrusion length is zero.
        if not salty.any():
            front_x[t] = x_mouth
            intrusion_length[t] = 0.0
            continue

        # Entire transect is salty: salt reaches the upstream boundary.
        if salty.all():
            front_x[t] = x_upstream
            intrusion_length[t] = x_mouth - x_upstream
            continue

        # Search from the mouth upstream, without reversing arrays.
        #
        # We look for the first sustained fresh segment encountered while
        # moving upstream from the mouth. The front is then the first salty
        # point immediately downstream of that sustained fresh segment.
        fresh_count = 0
        seen_salty_downstream = False
        front_idx = None

        for i in range(n_reg - 1, -1, -1):  # mouth -> upstream
            if salty[i]:
                seen_salty_downstream = True
                fresh_count = 0
            else:
                if seen_salty_downstream:
                    fresh_count += 1

                    if fresh_count >= min_consecutive:
                        # Sustained fresh segment found.
                        # The first salty point downstream should be:
                        candidate = i + min_consecutive

                        if candidate < n_reg and salty[candidate]:
                            front_idx = candidate
                        else:
                            # Robust fallback: find nearest salty point downstream.
                            downstream_salty = np.where(salty[i + 1:])[0]
                            if downstream_salty.size > 0:
                                front_idx = i + 1 + downstream_salty[0]

                        break

        if front_idx is None:
            # Fallback: use the most upstream salty point.
            front_idx = np.where(salty)[0][0]

        # Optional sub-grid interpolation of the threshold crossing.
        #
        # In the ideal case:
        #   front_idx - 1 is fresh
        #   front_idx     is salty
        #
        # Then interpolate between them to estimate where salinity = threshold.
        if front_idx > 0 and profile[front_idx - 1] < threshold <= profile[front_idx]:
            x0 = x_reg[front_idx - 1]
            x1 = x_reg[front_idx]
            s0 = profile[front_idx - 1]
            s1 = profile[front_idx]

            if s1 != s0:
                xf = x0 + (threshold - s0) * (x1 - x0) / (s1 - s0)
            else:
                xf = x_reg[front_idx]
        else:
            xf = x_reg[front_idx]

        front_x[t] = xf
        intrusion_length[t] = x_mouth - xf

    # ---------------------------------------------------------
    # 5. Optional temporal smoothing
    # ---------------------------------------------------------
    if temporal_smooth_hours is not None and dt_hours is not None:
        win_t = max(1, int(round(temporal_smooth_hours / dt_hours)))

        intrusion_length = uniform_filter1d(
            intrusion_length,
            size=win_t,
            mode="nearest",
        )

        # Smooth front_x consistently if requested.
        front_x = x_mouth - intrusion_length

    return intrusion_length


def plot_intrusion_time_series(x, salinity, threshold=0.5):
    """
    Plot time series of salt intrusion length along a transect.
    x: 1D array of distance along transect (e.g., river miles)
    salinity: 2D array of salinity (nt, npoints)
    threshold: salinity threshold to define intrusion
    """
    intrusion_length = cal_intrusion_length_from_mouth(x, salinity, threshold)

    plt.figure(figsize=(10, 4))
    plt.plot(intrusion_length, marker='o')
    plt.xlabel('Time Index')
    plt.ylabel('Salt Intrusion Length (distance units)')
    plt.title(f'Salt Intrusion Length Time Series (Threshold={threshold})')
    plt.grid()
    plt.show()


def compare_intrusion_lengths(threshold=0.5, ramp_up_t_records=7*24, reverse=False, use_miles=False):
    """
    Compare salt intrusion lengths from two different salinity datasets.
    x: 1D array of distance along transect (e.g., river miles), must be from downstream to upstream
    salinity1, salinity2: 2D arrays of salinity (nt, npoints)
    threshold: salinity threshold to define intrusion
    ramp_up_t_records: number of initial time records to skip for ramp-up
    reverse: if True, reverse the x array to be from upstream to downstream before calculation
    use_miles: if True, convert intrusion length to miles; otherwise, use km
    """
    transect_name = 'FreshwaterBayouCanal' # 'mississippi'  # 'FreshwaterBayouCanal'
    results_list = {
        'coarse mesh': f'/sciclone/schism10/feiye/STOFS3D-v8/O29i1/salinity.transect.{transect_name}.nc',
        'refined mesh': f'/sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.transect.{transect_name}.nc',
    }
    start_time = [pd.Timestamp('2024-03-05', tz='UTC'), pd.Timestamp('2024-03-05', tz='UTC')]
    time_stamps = []
    for result in results_list.values():
        ds = xr.open_dataset(result)
        time_values = ds['time'].values  # (nt,)
        time_datetimes = start_time[0] + pd.to_timedelta(time_values, unit='days')
        time_stamps.append(time_datetimes)

    intrusion_data_list = []
    for time_stamp, (label, filepath) in zip(time_stamps, results_list.items()):
        ds = xr.open_dataset(filepath)
        time_stamp = time_stamp[ramp_up_t_records:]
        salinity = ds['salinity'].values  # (nt, npoints)
        salinity = salinity[ramp_up_t_records:, :, 0]  # bottom
        x = compute_along_transect_distance(ds['lon'].values, ds['lat'].values, reverse=reverse)  # in miles
        x *= 1609.34  # miles to meters

        intrusion_length = cal_intrusion_length_from_mouth(x, salinity, threshold)  # meters from mouth
        if use_miles:
            intrusion_length /= 1609.34  # convert to miles
        else:
            intrusion_length /= 1000.0  # convert to km
        intrusion_data_list.append((intrusion_length))

        ylim = [np.min(intrusion_length) - 8, np.max(intrusion_length) + 6]
        plt.plot(time_stamp, intrusion_length, label=label)
        plt.ylim(ylim[0], ylim[1])
    
    mean_diff = np.nanmean(np.abs(
        intrusion_data_list[0] - intrusion_data_list[1]
    ))
    
    # plot
    ylable_str = 'Salt intrusion length (miles)' if use_miles else 'Salt intrusion length (km)'
    plt.ylabel(ylable_str)
    plt.title(f'Comparison of Salt Intrusion Lengths (Threshold={threshold}); Mean Abs Diff={mean_diff:.2f} miles)')
    plt.legend()
    plt.grid()
    plt.show()


def compute_along_transect_distance(lon, lat, reverse=False):
    """
    Compute cumulative distance along a transect given lon/lat coordinates.
    Returns distance in miles.
    """
    # project lon/lat to esri:102008
    import geopandas as gpd
    gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon, lat), crs='EPSG:4326')
    gdf = gdf.to_crs('ESRI:102008')
    river_x, river_y = gdf.geometry.x.values, gdf.geometry.y.values
    river_dist = np.cumsum(np.sqrt(np.diff(river_x, prepend=river_x[0])**2 + np.diff(river_y, prepend=river_y[0])**2))
    river_dist -= river_dist[0]  # start from 0
    river_dist_miles = river_dist / 1609.34  # meters to miles

    if reverse:
        river_dist_miles = river_dist_miles[-1] - river_dist_miles
    return river_dist_miles
    

def plot_model_transect2D(
    processed_schism_outputs='/sciclone/schism10/feiye/STOFS3D-v8/O29j/salinity.transect.mississippi.nc',
    schism_start_time=None,  # todo: put this into processed netcdf
    var_str='salinity', plot_time=pd.Timestamp('2023-08-31', tz='UTC'),
    ax=None, clim=None, draw_mesh=False, reverse_x=False,
    draw_color_bar=True, river_mile_coor=True, show_plot=True, use_miles=False
):
    """
    Plot a 2D transect from SCHISM model output along a river,
    with river miles added and off-limit miles filled in.
    """
    if schism_start_time is None:
        raise ValueError("schism_start_time must be provided.")
    results_ds = xr.open_dataset(processed_schism_outputs, engine='netcdf4')
    run_name = Path(processed_schism_outputs).parent.name

    # find nearest time index to plot_time
    time_values = results_ds['time'].values  # (nt,)
    time_datetimes = schism_start_time + pd.to_timedelta(time_values, unit='days')
    time_deltas = np.abs(time_datetimes - plot_time)
    t_idx = np.argmin(time_deltas)
    print(f'Plotting transect at the nearest time: {time_datetimes[t_idx]} (index {t_idx})')

    var_name_dict = {
        'salinity': 'Salinity (psu)',
        'temperature': 'Temperature (°C)',
        'diffusivity': 'Diffusivity (m²/s)',
        'viscosity': 'Viscosity (kg/(m·s))',
    }

    river_lon = results_ds['lon'].values  # (npoints,)
    river_lat = results_ds['lat'].values  # (npoints,)
    z = results_ds['zCoordinates'].values  # (nt, npoints, nz)
    var = results_ds[var_str].values   # (nt, npoints, nz)

    valid_col = z[0, :, 0] > -9000  # find columns with NaN
    river_lon = river_lon[valid_col]
    river_lat = river_lat[valid_col]
    z = z[:, valid_col, :]
    var = var[:, valid_col, :]
    
    # compute river miles along the transect
    if river_mile_coor:
        river_miles_ref = RiverMilesRef.from_file(mile_max=300, mile_min=-10)
        miles = np.array([river_miles_ref.to_mile(lon, lat) for lon, lat in zip(river_lon, river_lat)], dtype=np.float32)

    if river_mile_coor:
        # compute along transect distance to fill in off-limit miles
        river_dist_miles = compute_along_transect_distance(river_lon, river_lat)
        # fill off-limit miles with river distance in miles
        upstream_invalid_idx = miles >= 300 - 1 
        downstream_invalid_idx = miles <= 0
    
        if np.any(upstream_invalid_idx):
            first_valid_idx = np.where(upstream_invalid_idx)[0][-1] + 1
            upstream_invalid_idx = np.arange(0, first_valid_idx)  # also fill any earlier valid points mixed with invalids
            fill_miles_increments = river_dist_miles[first_valid_idx] - river_dist_miles[upstream_invalid_idx]
            miles[upstream_invalid_idx] = miles[first_valid_idx] + fill_miles_increments

        if np.any(downstream_invalid_idx):
            first_valid_idx = np.where(downstream_invalid_idx)[0][0] - 1
            downstream_invalid_idx = np.arange(first_valid_idx + 1, len(miles))  # also fill any later valid points mixed with invalids
            fill_miles_decrements = river_dist_miles[downstream_invalid_idx] - river_dist_miles[first_valid_idx]
            miles[downstream_invalid_idx] = miles[first_valid_idx] - fill_miles_decrements
        
        plot_transect_filled(
            miles, -z, var, var_name_dict[var_str], t_idx=t_idx,
            ax=ax, cmap='jet', clim=clim, draw_mesh=draw_mesh, draw_color_bar=draw_color_bar
        )
    
    else:
        # reverse river_dist_miles if necessary, start from downstream
        river_dist_miles = compute_along_transect_distance(river_lon, river_lat, reverse=False)
        river_dist_miles = river_dist_miles[-1] - river_dist_miles  # from downstream to upstream
        if use_miles:
            river_dist = river_dist_miles  # already in miles
        else:
            river_dist = river_dist_miles * 1.60934  # miles to km

        if var_str == 'temperature':
            clim = (15.5, 22.5)
        elif var_str == 'salinity':
            clim = (0, 4)
        ax = plot_transect_filled(
            river_dist, -z, var, var_name_dict[var_str], ax=ax, t_idx=t_idx, cmap='jet',
            draw_mesh=draw_mesh, reverse_x=reverse_x, draw_color_bar=draw_color_bar,
            clim=clim,
            # plot_args={'xlim': (0, 45*1.60934)}
        )

        ax.set_ylabel("z (m)")
        if not use_miles:
            ax.set_xlabel("along-transect distance (km)")  # override default label for x-axis
        plt.title(f"Model {var_str.capitalize()} Transect — {run_name} at {time_datetimes[t_idx].date()}")

    if show_plot:
        plt.show()
    pass

    return time_datetimes[t_idx], t_idx


def diff_transect_2D(var_str='salinity', x_range=None, schism_start_time=None, plot_time=None):
    """
    Plot difference between two 2D transects from SCHISM model outputs.

    Input:
    ------
    var_str: str
        Variable to plot difference for ('salinity' or 'temperature').
    x_range: tuple or None
        x-axis range (in miles) to plot. If None, plot full range.
    schism_start_time: pd.Timestamp or None
    plot_time: pd.Timestamp or None
        Time to plot. If None, plot at t_idx=0.
    """
    transect_name = 'FreshwaterBayouCanal'  # 'mississippi'  # 'FreshwaterBayouCanal'
    ds1 = xr.open_dataset(f'/sciclone/schism10/feiye/STOFS3D-v8/O19i1/salinity.transect.{transect_name}.nc', engine='netcdf4')
    ds2 = xr.open_dataset(f'/sciclone/schism10/feiye/STOFS3D-v8/O29i1/salinity.transect.{transect_name}.nc', engine='netcdf4')

    var_name_dict = {
        'salinity': 'Salinity Difference (psu)',
        'temperature': 'Temperature Difference (°C)',
    }

    if plot_time is not None and schism_start_time is not None:
        time_values = ds1['time'].values  # (nt,)
        time_datetimes = schism_start_time + pd.to_timedelta(time_values, unit='days')
        time_deltas = np.abs(time_datetimes - plot_time)
        t_idx = np.argmin(time_deltas)
        print(f'Plotting transect difference at the nearest time: {time_datetimes[t_idx]} (index {t_idx})')
    elif plot_time is None and schism_start_time is None:
        t_idx = 0
    else:
        raise ValueError("Both schism_start_time and plot_time must be provided together.")

    river_lon = ds1['lon'].values  # (npoints,)
    river_lat = ds1['lat'].values  # (npoints,)
    if not np.allclose(river_lon, ds2['lon'].values) or not np.allclose(river_lat, ds2['lat'].values):
        raise ValueError("The two datasets have different transect coordinates.")

    z1 = ds1['zCoordinates'].values  # (nt, npoints, nz)
    var1 = ds1[var_str].values   # (nt, npoints, nz)
    z2 = ds2['zCoordinates'].values  # (nt, npoints, nz)
    var2 = ds2[var_str].values   # (nt, npoints, nz)

    valid_col = z1[0, :, 0] > -9000  # find columns with NaN
    river_lon = river_lon[valid_col]
    river_lat = river_lat[valid_col]
    z1 = z1[:, valid_col, :]
    var1 = var1[:, valid_col, :]
    z2 = z2[:, valid_col, :]
    var2 = var2[:, valid_col, :]

    var_diff = var1 - var2

    river_dist_miles = compute_along_transect_distance(river_lon, river_lat, reverse=False)
    river_dist_miles = river_dist_miles[-1] - river_dist_miles
    if x_range is not None:
        mask = (river_dist_miles >= x_range[0]) & (river_dist_miles <= x_range[1])
        river_dist_miles = river_dist_miles[mask]
        z1 = z1[:, mask, :]
        var_diff = var_diff[:, mask, :]

    rmsd = np.sqrt(np.nanmean(var_diff[t_idx]**2))

    plot_transect_filled(
        river_dist_miles, -z1, var_diff, var_name_dict[var_str], t_idx=t_idx,
        cmap='bwr', clim=(-2*rmsd, 2*rmsd), draw_mesh=True, reverse_x=False,
        draw_color_bar=True,
    )
    plt.title(
        f"Difference in {var_str.capitalize()} Transect at t_idx={t_idx}; "
        f"RMSD={rmsd:.2f}"
    )
    plt.show()

if __name__ == "__main__":
    # set global plot parameters
    plt.rcParams.update({'font.size': 14})

    # compare_intrusion_lengths(threshold=4, reverse=False, use_miles=False)

    # diff_transect_2D(
    #     var_str='salinity',
    #     # x_range=(0, 45),
    #     schism_start_time=pd.Timestamp('2024-03-05', tz='UTC'),
    #     plot_time=pd.Timestamp('2024-04-09', tz='UTC')
    # )

    plot_model_transect2D(
        processed_schism_outputs='/sciclone/schism10/feiye/STOFS3D-v8/O29i1/salinity.transect.FreshwaterBayouCanal.nc',
        schism_start_time=pd.Timestamp('2024-03-05', tz='UTC'),
        var_str='salinity', plot_time=pd.Timestamp('2024-04-09 00:00:00', tz='UTC'),
        river_mile_coor=False, draw_mesh=False, show_plot=True, use_miles=False
    )

    pass
