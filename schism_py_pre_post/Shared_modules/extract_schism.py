import shutil
import os
import copy
import subprocess


def extract_schism(bpfiles=None, var=None, rundir=None, i_comb=None, stacks=None, i_overwrite=True, i_all_level=False, ver=9):
    '''
    extract schism results based on one or more bpfiles
    ver: output version (e.g., ver 10 needs read_output10*)
    '''
    runid = os.path.basename(os.path.normpath(rundir))

    if i_all_level:
        script = f'/sciclone/home/feiye/bin/read_output{ver}_xy_all_level.viz'
    else:
        script = f'/sciclone/home/feiye/bin/read_output{ver}_xyz.viz'

    out_filenames = []

    if not os.path.exists(rundir):
        raise Exception("Schism dir not exist: " + rundir)
    if not os.path.exists(rundir + '/outputs/'):
        raise Exception("Schism outputs dir not exist: " + rundir)
    if not os.path.exists(rundir + '/extract_in/'):
        os.mkdir(rundir + "/extract_in/")
    if os.path.islink(rundir + '/PostP'):
        os.remove(rundir + '/PostP')
    if not os.path.exists(rundir + '/PostP/'):
        os.mkdir(rundir + "/PostP/")
    if not os.path.exists(rundir + '/outputs/vgrid.in'):
        shutil.copy2(rundir + '/vgrid.in', rundir + '/outputs/')

    for bpfile in bpfiles:
        bpname = os.path.basename(os.path.normpath(bpfile))
        bpname = bpname.split('.')[0]
        shutil.copy2(bpfile, rundir + '/outputs/station.bp')

        # make read_output*.in
        read_in = rundir + '/extract_in/' + var + '.in'
        with open(read_in, 'w') as fout:
            if ver==9:
                if i_comb:
                    fout.write("1\n2\n10\n")
                else:
                    fout.write("0\n1\n2.e9\n")
            fout.write("1\n1\n")
            fout.write(var + '\n')
            fout.write("1\n")
            fout.write(str(stacks[0]) + ' ' + str(stacks[1]) + '\n')
            fout.write("1\n")

        # call fortran script
        if i_all_level:
            out_filename = rundir + '/PostP/' + var + '.all_level.' + bpname + '.' + runid
        else:
            out_filename = rundir + '/PostP/' + var + '.dat.' + bpname + '.' + runid

        if var == 'hvel':
            out_filenames.append(out_filename.replace('hvel', 'u'))
            out_filenames.append(out_filename.replace('hvel', 'v'))
        else:
            out_filenames.append(out_filename)

        for this_file in out_filenames:
            if not os.path.exists(this_file) or i_overwrite:
                print(">>>> extracting schism results ...")
                p = subprocess.Popen(script, stdin=open(read_in), cwd=(rundir + '/outputs/'))
                p.wait()
                # input("Press Enter to continue...")
                if var == 'hvel':
                    os.system(f'mv {rundir}/outputs/fort.18 {out_filename.replace("hvel", "u")}')
                    os.system(f'mv {rundir}/outputs/fort.19 {out_filename.replace("hvel", "v")}')
                else:
                    os.system(f'mv {rundir}/outputs/fort.18 {out_filename}')
            else:
                print(">>>> schism results already extracted, skipping ...")

    return out_filenames

if __name__ == "__main__":
    '''Samples'''

    # ------------------------- time series --------------------------- 
    # rundir = '/sciclone/data10/feiye/vims20/work/ChesBay/RUN200p/'
    # bpname = 'center'
    # outfilenames = extract_schism(
    #     bpfiles=[f'/sciclone/data10/feiye/vims20/work/ChesBay/BPfiles/{bpname}.bp'],
    #     var='temperature',
    #     rundir=rundir,
    #     i_comb=True,
    #     stacks=[1, 255],
    #     i_all_level=False,
    #     ver=10
    # )
    # shutil.copy(outfilenames[0], f'{rundir}/temp.dat.more.{bpname}')
    # pass

    # # ------------------------- ChesBay transect (synoptic_transect3.m) --------------------------- 
    # bpfiles = ['/sciclone/data10/feiye/vims20/work/ChesBay/BPfiles/side_w.bp']
    # 
    # runid = 'RUN109x'
    # vars = ['salt', 'zcor']
    # stacks = [1, 86]
    # ver=9
    # 
    # runid = 'RUN200p'
    # vars = ['salinity', 'zCoordinates']
    # stacks = [1, 255]
    # ver=10
    # 
    # rundir = f'/sciclone/data10/feiye/vims20/work/ChesBay/{runid}/'
    # for var in vars:
    #     outfilenames = extract_schism(
    #         bpfiles=bpfiles,
    #         var=var,
    #         rundir=rundir,
    #         i_comb=True,
    #         stacks=stacks,
    #         i_all_level=True,
    #         ver=ver
    #     )
    # pass

    # ------------------------- ChesBay hvel --------------------------- 
    # runid = 'RUN200i'
    # rundir = '/sciclone/data10/feiye/vims20/work/ChesBay/' + runid
    # bpdir = '/sciclone/data10/feiye/vims20/work/ChesBay/BPfiles/'
    # bpnames = ['cb0102','cb0301','cb0601','cb0701','cb0901','cb1001','cb1101','cb1201','cb0402']
    # for bpname in bpnames:
    #     for var1, var2 in [['horizontalVelX', 'u'], ['horizontalVelY', 'v']]:
    #         outfilenames = extract_schism(
    #             bpfiles=[f'{bpdir}/{bpname}.bp'],
    #             var=var1,
    #             rundir=rundir,
    #             i_comb=True,
    #             stacks=[1, 80],
    #             i_all_level=False,
    #             ver=10
    #         )
    #         shutil.copy(outfilenames[0], f'{rundir}/{var2}.dat.{bpname}.{runid}')
    #     pass

    # extract_schism(
    #     bpfiles=['/sciclone/schism10/feiye/ICOGS/BPfiles/missi.bp'],
    #     var='elev',
    #     rundir='/sciclone/schism10/feiye/ICOGS/RUN22b/',
    #     i_comb=False,
    #     stacks=[1, 24],
    #     i_all_level=False
    # )


    rundir = '/sciclone/schism10/feiye/STOFS3D-v7/Outputs/O15/'
    bpname = "Missi_Ida2_moved"
    outfilenames = extract_schism(
        bpfiles=[f'/sciclone/schism10/feiye/ICOGS/BPfiles/{bpname}.bp'],
        var='elevation',
        rundir=rundir,
        i_comb=True,
        stacks=[1, 55],
        i_all_level=False,
        ver=10
    )
    print(f'extraction done: {outfilenames}')
