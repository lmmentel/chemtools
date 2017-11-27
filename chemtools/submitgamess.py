# -*- coding: utf-8 -*-

#The MIT License (MIT)
#
#Copyright (c) 2014 Lukasz Mentel
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import argparse
import os
import socket
import subprocess


def main():
    '''
    Script for submitting gamessus jobs to the queue.
    '''

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("input",
                        help='gamessus input file to be executed')
    parser.add_argument("-d",
                        "--dryrun",
                        action="store_true",
                        help="write the pbs script but don't submit to the queue, default-False")
    parser.add_argument("-e",
                        "--extrafiles",
                        nargs="+",
                        default=[],
                        help="additional files that need to be copied to the node/scratch")
    parser.add_argument("-g",
                        "--gmsver",
                        default="00",
                        help="version of gamess executable, default=00")
    group.add_argument("-n",
                       "--nodes",
                       default="1",
                       help="number of nodes, default=1")
    group.add_argument("-H",
                       "--HOST",
                       default="",
                       help="destination hostname, default=''")
    parser.add_argument("-p",
                        "--ppn",
                        default="1",
                        help="number of processes per node, default=1")
    parser.add_argument("-s",
                        "--usescratch",
                        action="store_true",
                        help="use node scratch directory (/scratch/$USER)",
                        default=False)
    parser.add_argument("-q",
                        "--queue",
                        default="default",
                        help="destination queue, default=default")
    parser.add_argument("-t",
                        "--walltime",
                        default="120:00:00",
                        help="walltime in the format HH:MM:SS")
    args = vars(parser.parse_args())
    submit_pbs(args)


def set_defaults(args):

    args['workdir'] = os.getcwd()
    args['home'] = os.getenv("HOME")
    if socket.gethostname() == "boron.chem.umk.pl":
        args['scratch'] = os.path.join('/scratch', os.getenv('USER'))
        args['rungms'] = '/share/apps/gamess-may2013/rungms'
    elif socket.gethostname() in ["login1.lisa.surfsara.nl",
                                  "login2.lisa.surfsara.nl"]:
        args['scratch'] = os.getenv("TMPDIR")
        args['rungms'] = 'rungms'
    args['local_scr'] = os.path.join(os.getenv('HOME'), 'scratch')
    args['jobname'] = os.path.splitext(args["input"])[0]
    args['outfile'] = args['jobname'] + ".log"
    args['errfile'] = args['jobname'] + ".err"
    args['datfile'] = args['jobname'] + ".dat"
    args['script_name'] = "run." + args['jobname']
    return args


def remove_dat(path, datfile):
    '''Remove the dat file from the ASCII scratch.'''
    if os.path.exists(os.path.join(path, datfile)):
        os.remove(os.path.join(path, datfile))


def submit_pbs(args):
    '''
    Write the run script for PBS and submit it to the queue.
    '''

    args = set_defaults(args)
    remove_dat(args["local_scr"], args["datfile"])

    with open(args['script_name'], 'w') as script:
        script.write("#PBS -S /bin/bash\n")
        if args['HOST'] != "":
            script.write("#PBS -l nodes={0}:ppn={1}\n".format(args['HOST'], args['ppn']))
        else:
            script.write("#PBS -l nodes={0}:ppn={1}\n".format(args['nodes'], args['ppn']))
        if "mem" in args.keys():
            script.write("#PBS -l mem={0}\n".format(args["mem"]))
        script.write("#PBS -l walltime={0}\n\n".format(args['walltime']))
        #script.write('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/lib\n\n')
        script.write("#PBS -N {}\n".format(args["jobname"]))
        if args["queue"] != "default":
            script.write("#PBS -q {}\n".format(args["queue"]))
        script.write("cd $PBS_O_WORKDIR\n")
        if args['usescratch']:
            wrkdir = os.path.join(args['scratch'], args['jobname'])
            script.write("mkdir -p {}\n".format(wrkdir))
            files = args['jobname']
            if args['extrafiles']:
                files += ' ' + ' '.join(args['extrafiles'])
            script.write('cp -t {0} {1}\n'.format(wrkdir, files))
            script.write('cd {0}\n'.format(wrkdir))
        script.write("\n{0:<s} {1} {2:<s} {3:<s} \n".format(args['rungms'],
                     args['jobname'], args['gmsver'], args['ppn']))

    # submit the job to the queue if requested
    if args['dryrun']:
        print("Created job script: {0}\n NOT submitting to the queue\nbye...".format(args['script_name']))
    else:
        print("Created job script: {0}\nsubmitting to the queue".format(args['script_name']))
        output = subprocess.check_output(["qsub", args['script_name']])
        pid = output.split(".")[0]
        return args['jobname'] + ".o" + pid


def submit_slurm(args):
    '''
    Write the run script for SLURM and submit it to the queue.
    '''

    with open(args['script_name'], 'w') as script:

        script.write("#!/bin/bash\n")
        script.write("#SBATCH -t {0:<s} \n".format(args["walltime"]))
        script.write("#SBATCH -N {0:<s} \n".format(args["nodes"]))
        script.write("#SBATCH --ntasks-per-node={0:<s} \n".format(args["ppn"]))
        script.write("#SBATCH -C {0:<s} \n".format(args["nodeType"]))
        script.write("#SBATCH -o {0:>s}.%J \n".format(args['outfile']))
        script.write("#SBATCH -e {0:>s}.%J \n".format(args['errfile']))

        if int(args["walltime"].split(':')[0])*60 + int(args["walltime"].split(':')[1]) > 60 and args["nodeType"] == "thin":
            script.write("#SBATCH -p normal \n\n")
        elif int(args["walltime"].split(':')[0])*60 + int(args["walltime"].split(':')[1]) < 60 and args["nodeType"] == "thin":
            script.write("#SBATCH -p short \n\n")
        elif args["nodeType"] == "fat":
            script.write("#SBATCH -p fat \n\n")

        if args["mail"] != "":
            script.write("#SBATCH --mail-type=ALL --mail-user={0:<s}\n\n".format(args["mail"]))
        script.write("mkdir -p $TMPDIR/{0}\n".format(args['jobname']))
        script.write("cp {0} $TMPDIR/{1} \n".format(args["input"], args['jobname']))
        script.write("cd $TMPDIR/{0}\n\n".format(args['jobname']))

        script.write("module load gamess-us \n")
        script.write("{0} {1} \n".format(args['rungms'], args['jobname']))
        #script.write("cp {0} {1}\n".format(output, workdir))

    # submit the job to the queue if requested
    if args["dryrun"]:
        print("Created job script: {0}\n NOT submitting to the queue\nbye...".format(args['script_name']))
    else:
        print("Created job script: {0}\nsubmitting to the queue".format(args['script_name']))
        sublog = open(args['jobname'] + ".sublog", 'w')
        proc = subprocess.Popen(["sbatch", args['script_name']], stdout=sublog, stderr=sublog)
        sublog.close()


def submit_ll(args):
    '''
    Write the run script for LoadLeveller and submit it to the queue.
    '''

    with open(args['script_name'], 'w') as script:

        if args['ppn'] > 1 or args['nodes'] > 1:
            job_type = "parallel"
            script.write("# @ node = {0}\n".format(args['nodes']))
            script.write("# @ tasks_per_node = {0}\n".format(args['ppn']))
            script.write("# @ network.MPI_LAPI = sn_all,not_shared,US \n")
        else:
            job_type = "serial"
            script.write("# @ node_usage = shared \n")

        script.write("#\n# @ notification = never\n#\n# @ input = /dev/null\n")
        script.write("# @ output = {0}.$(jobid)\n".format(args['outfile']))
        script.write("# @ error = {0}.$(jobid)\n".format(args['errfile']))
        script.write("# @ wall_clock_limit = {0}\n#\n".format(args['walltime']))
        script.write("# @ job_type = {0}\n#\n".format(job_type))
        script.write("# @ queue \n#\n")
        script.write("mkdir -p {0} \n".format(args['scratch']))
        script.write("cp {0} {1} \n".format(
            os.path.join(args['workdir'], args['jobname']),
            os.path.join(args['scratch'], args['jobname'])))
        script.write("cd {0} \n".format(args['scratch']))
        script.write("module load gamess-us/May2012 \n")
        script.write("{0} {1} \n".format(args['rungms'], args['jobname']))

    # submit the job to the queue if requested
    if args["dryrun"]:
        print("Created job script: {0}\n NOT submitting to the queue\nbye...".format(args['script_name']))
    else:
        print("Created job script: {0}\nsubmitting to the queue".format(args['script_name']))
        sublog = open(args['jobname'] + ".sublog", 'w')
        proc = subprocess.Popen(["llsubmit", args['script_name']], stdout=sublog, stderr=sublog)
        sublog.close()


if __name__ == "__main__":
    main()
