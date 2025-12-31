#!/usr/bin/env python3

"""
Make a new run directory by linking input files from another
"""


import os
import argparse


def make_new_run(old_run_dir, new_run_dir, scr):
    '''
    Make a new run by:
    1) link input files from old_run_dir to new_run_dir;
    2) make a new outputs directory in on scratch and link to new_run_dir.
    3) copy some of the input files to the new run directory
    '''

    os.makedirs(new_run_dir)
    os.chdir(new_run_dir)

    # link input files
    os.system(f'ln -s ../{old_run_dir}/* .')

    # remove outputs
    os.system('rm -rf *out* *dat*')

    # make a new outputs directory in on scratch and link to new_run_dir
    if scr is None:
        os.makedirs('outputs/', exist_ok=True)
    else:
        os.makedirs(f'{scr}/{new_run_dir}/outputs/', exist_ok=True)
        os.system(f'ln -s {scr}/{new_run_dir}/outputs/ .')

    # copy some of the input files to the new run directory
    os.system('rm param.nml run_test')
    os.system(f'cp ../{old_run_dir}/param.nml .')
    os.system(f'cp ../{old_run_dir}/run_test .')


def handle_cmd_args():
    '''
    Handle command line arguments
    '''
    parser = argparse.ArgumentParser(
        description='Make a new run directory by linking input files from another')
    parser.add_argument('old_run_dir', type=str, help='old run directory')
    parser.add_argument('new_run_dir', type=str, help='new run directory')
    parser.add_argument('scr', type=str, nargs="?", help='scratch directory', default=None)
    args = parser.parse_args()

    print(f'old_run_dir: {args.old_run_dir}')
    print(f'new_run_dir: {args.new_run_dir}')
    print(f'scr: {args.scr}')

    return args


def main():
    '''
    Main function
    '''
    args = handle_cmd_args()
    make_new_run(args.old_run_dir, args.new_run_dir, args.scr)


if __name__ == '__main__':
    main()

