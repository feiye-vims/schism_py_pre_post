'''
Utility functions for general purposes.
'''


from datetime import datetime
import subprocess
import numpy as np
from tqdm import tqdm


def b_in_a(a=None, b=None):
    '''
    Given two arrays A and B, return the indices of A's elements in B.
    '''

    if a is None and b is None:
        print('Demonstration of b_in_a')
        a = np.array([3, 5, 7, 1, 9, 8, 6, 6])
        b = np.array([3, 1, 5, 8, 6])
        print(f'Array a: {a}')
        print(f'Array b: {b}')
    else:
        a = np.array(a)
        b = np.array(b)

    index = np.argsort(a)
    sorted_a = a[index]
    sorted_index = np.searchsorted(sorted_a, b)

    bindex = np.take(index, sorted_index, mode="raise")
    mask = a[bindex] != b

    result = np.ma.array(bindex, mask=mask)
    return result


def parse_date(date_string):
    '''
    Parses a date string and returns a datetime object.
    The order of the formats in the list matters:
    the function will stop at the first format that matches the input string.
    '''

    formats = [
        "%Y-%m-%d",           # e.g. 2023-06-13
        "%Y-%m-%d %H:%M:%S",  # e.g. 2023-06-13 00:00:00
        "%Y-%m-%d %H:%M",     # e.g. 2023-06-13 00:00
        "%Y-%m-%d %H",        # e.g. 2023-06-13 00
        "%Y-%m-%dT%H:%M:%S",  # e.g. 2023-06-13 00:00:00
        "%Y-%m-%dT%H:%M",     # e.g. 2023-06-13 00:00
        "%Y-%m-%dT%H",        # e.g. 2023-06-13 00
        "%m/%d/%Y",           # e.g. 06/13/2023
        "%d/%m/%Y",           # e.g. 13/06/2023
        "%B %d, %Y",          # e.g. June 13, 2023
        "%d %B, %Y",          # e.g. 13 June, 2023
        "%m-%d-%Y",           # e.g. 06-13-2023
        "%d-%m-%Y",           # e.g. 13-06-2023
        "%Y/%m/%d",           # e.g. 2023/06/13
        "%m.%d.%Y",           # e.g. 06.13.2023
        "%d.%m.%Y",           # e.g. 13.06.2023
        "%Y.%m.%d",           # e.g. 2023.06.13
        "%b %d, %Y",          # e.g. Jun 13, 2023
        "%d %b, %Y",          # e.g. 13 Jun, 2023
        "%A %B %d, %Y",       # e.g. Tuesday June 13, 2023
        "%d %B %A, %Y",       # e.g. 13 June Tuesday, 2023
        "%d %b %a, %Y",       # e.g. 13 Jun Tue, 2023
        "%d-%b-%y",           # e.g. 13-Jun-23
    ]

    for fmt in formats:
        try:
            return datetime.strptime(date_string, fmt), fmt
        except ValueError:
            pass
    print("--------------Invalid date format!----------------")
    return None, None


def my_mpi_idx(n_tasks, size, rank, distribution_type='cyclic'):
    '''
    Distribute N tasks to {size} ranks.
    The return value is a bool vector of the shape (n_tasks, ),
    with True indices indicating tasks for the current rank.
    '''

    # initialize a bool array of size n_tasks
    # for the current rank to decide which tasks to handle
    # (True for the tasks to handle, False for the tasks to skip)
    my_idx = np.zeros((n_tasks, ), dtype=bool)

    if n_tasks <= size:  # trivial case: more ranks than tasks
        if rank < n_tasks:
            my_idx[rank] = True
    else:
        if distribution_type == 'cyclic':
            for i in range(n_tasks):
                if (i % size) == rank:
                    my_idx[i] = True
        elif distribution_type == 'block':
            my_idx = np.zeros((n_tasks, ), dtype=bool)
            n_per_rank, _ = divmod(n_tasks, size)
            n_per_rank = n_per_rank + 1
            my_idx[rank*n_per_rank:min((rank+1)*n_per_rank, n_tasks)] = True

    return my_idx


vdatum_preset = {
    'navd88_to_xgeoid20b':
        'ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height',
    'xgeoid20b_to_navd88':
        'ihorz:igs14:geo:deg ivert:xgeoid20b:m:height ohorz:NAD83_2011 overt:navd88:m:height',
}


def vdatum_wrapper_pointwise(x, y, z, conversion_para='', print_info=''):
    '''
    Wrapper function to convert points using VDatum
    The default conversion is from NAVD88 to XGeoid20b based on height

    Sample command for vdatum:
    {path_to_java}/java -jar vdatum.jar
    ihorz:NAD27:geo:deg ivert:DTL:US_ft:height ohorz:NAD83_2011:geo:deg overt:NAVD88:m:height
    -deg2dms -pt:-97.30965,26.3897528,3.545 region:4
    '''
    vdatum_folder = '/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum/'

    # z_convention = 'height'  # "sounding": positive downwards; "height": positive upwards

    # default conversion parameters
    if conversion_para == '':
        conversion_para = vdatum_preset['navd88_to_xgeoid20b']

    z_converted = z.copy()
    regions = [4, 5]

    for i in tqdm(range(len(x)), desc=f"{print_info} Processing"):
        success = False
        for _ in regions:
            result = subprocess.run(
                f"java -jar {vdatum_folder}/vdatum.jar"
                "{conversion_para} -pt:{x[i]},{y[i]},{z[i]} region:{region}",
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

            if result.returncode == 0:
                z_converted[i] = float(result.stdout.decode().split()[417])
                print(f"{print_info}Point {i+1} ({x[i]}, {y[i]}, {z[i]})"
                      "successfully converted to {z_converted[i]}")
                success = True
                continue  # found the correct region, skip the rest

        if not success:
            print(f"{print_info}Point {i+1} ({x[i]}, {y[i]}, {z[i]}) failed to convert")
            # print("swapping default regions")
            # regions = [regions[1], regions[0]]

        # set toloreance to 1e-4, i.e., 0.1 mm
        # if not np.isclose(z[i], float(result.stdout.decode().split()[416]), atol=1e-4):
        #    raise Exception(f"Input z and output z do not match for point:"
        #                    "{i+1} ({x[i]}, {y[i]}, {z[i]})")

    return z_converted
