#!/usr/bin/env python
"""
This module includes functions related to date and time.
"""
from datetime import datetime, timedelta


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

def calculate_date_after_days(start_time_str, num_days, keep_original_fmt=True):
    '''
    Calculates the date after a given number of days from a given start date.
    By default, the output date will be in the same format as the input date.
    Otherwise, the output date will be in the format of "%Y-%m-%d %H:%M:%S".
    '''

    start_time, date_fmt = parse_date(start_time_str)
    end_time = start_time + timedelta(days=num_days)

    if keep_original_fmt:
        return end_time.strftime(date_fmt)
    else:
        return end_time.strftime("%Y-%m-%d %H:%M:%S")


if __name__ == "__main__":
    import sys

    # Ensure we have enough command line arguments
    if len(sys.argv) < 3:
        print("Usage: ./cal_time_from.py start_time_str number_of_days_float")
        sys.exit(1)

    string_arg = sys.argv[1]  # The first command line argument
    float_arg = sys.argv[2]   # The second command line argument

    # Attempt to convert the second argument to a float
    try:
        float_arg = float(float_arg)
    except ValueError:
        print("The second argument must be a float.")
        sys.exit(1)

    today = calculate_date_after_days(string_arg, float_arg)
    print(today)
