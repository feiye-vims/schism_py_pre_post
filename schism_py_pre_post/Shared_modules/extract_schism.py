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
        script = f'/sciclone/home10/feiye/bin/read_output{ver}_xy_all_level.viz'
    else:
        script = f'/sciclone/home10/feiye/bin/read_output{ver}_xyz.viz'

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
                    fout.write("1\n2\n1\n")
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


# -----------------------------------------------------------------------------------
# ----------------------extract schism results-------------------------
# -----------------------------------------------------------------------------------
rundir = '/sciclone/data10/feiye/vims20/work/ChesBay/RUN200i/'
outfilenames = extract_schism(
    bpfiles=['/sciclone/data10/feiye/vims20/work/ChesBay/BPfiles/center.bp'],
    var='salinity',
    rundir=rundir,
    i_comb=True,
    stacks=[1, 28],
    i_all_level=False,
    ver=10
)
shutil.copy(outfilenames[0], f'{rundir}/salt.dat.more.center')
pass

# extract_schism(
#     bpfiles=['/sciclone/schism10/feiye/ICOGS/BPfiles/missi.bp'],
#     var='elev',
#     rundir='/sciclone/schism10/feiye/ICOGS/RUN22b/',
#     i_comb=False,
#     stacks=[1, 24],
#     i_all_level=False
# )
