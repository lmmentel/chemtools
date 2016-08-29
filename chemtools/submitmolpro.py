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
    Script for submitting molpro job to the queue.
    '''

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("input",
                        help='molpro input file to be executed')
    parser.add_argument("-d",
                        "--dryrun",
                        action="store_true",
                        help="write the pbs script but don't submit to the queue, default-False")
    parser.add_argument("-e",
                        "--extrafiles",
                        nargs="+",
                        default=[],
                        help="additional files that need to be copied to the node/scratch")
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
                        help="walltime in the format HH:MM:SS,\
                        default=120:00:00")
    args = vars(parser.parse_args())
    submit_pbs(args)


def set_defaults(args):
    '''Set some useful default values and add them to args'''

    args['workdir'] = os.getcwd()
    args['home'] = os.getenv("HOME")
    if socket.gethostname() == "boron.chem.umk.pl":
        args['scratch'] = os.path.join('/scratch', os.getenv('USER'))
        args['molpro'] = '/share/apps/molprop_2012_1_Linux_x86_64_i8/bin/molpro'
    else:
        args['scratch'] = os.getenv("TMPDIR")
        args['molpro'] = 'molpro'
    args['local_scr'] = os.path.join(os.getenv('HOME'), 'scratch')
    args['jobname'] = os.path.splitext(args["input"])[0]
    args['outfile'] = args['jobname'] + ".o${PBS_JOBID}"
    args['errfile'] = args['jobname'] + ".e${PBS_JOBID}"
    args['script_name'] = "run." + args['jobname']
    return args


def submit_pbs(args):
    '''
    Write the run script for PBS and submit it to the queue.
    '''

    args = set_defaults(args)

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
        script.write("#PBS -o {}\n".format(args["outfile"]))
        script.write("#PBS -e {}\n".format(args["errfile"]))
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
        script.write("\n{0:<s} --nobackup -o {1} {2:s} \n".format(args['molpro'],
                     args['outfile'], args['jobname']))

    # submit the job to the queue if requested
    if args['dryrun']:
        print("Created job script: {0}\n NOT submitting to the queue\nbye...".format(args['script_name']))
    else:
        print("Created job script: {0}\nsubmitting to the queue".format(args['script_name']))
        output = subprocess.check_output(["qsub", args['script_name']])
        pid = output.split(".")[0]
        return args['jobname'] + ".o" + pid

if __name__ == "__main__":
    main()
