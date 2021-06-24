#!/usr/bin/env python
import os
import fileinput
import shutil
import json
import argparse
import sys
import subprocess

def reset_parameters():
    return {
        'FIELD_FILE'        : None,
        'FIELD_TOL'         : None,
        'H5_FIELD_FREQ'     : None,
        'H5_PARCEL_FREQ'    : None,
        'H5_BASENAME'       : None,
        'N_PER_CELL'        : None,
        'LAMBDA'            : None,
        'MERGE_TYPE'        : None,
        'VFRACTION'         : None,
        'VMAXFRACTION'      : None,
        'CORRECTION_ITERS'  : None,
        'GRADIENT_PREF'     : None,
        'MAX_COMPRESSION'   : None,
        'LIMIT'             : None,
        'DT'                : None,
        'ALPHA'             : None,
        'DT_MAX'            : None
    }


def load_json(fname):
    # 24 June 2021
    # https://stackoverflow.com/questions/20199126/reading-json-from-a-file
    with open(fname) as f:
        setup = json.load(f)
    return setup


def create_config_file(params, tmpl_dir, fname):
    # 24 June 2021
    # https://stackoverflow.com/questions/123198/how-can-a-file-be-copied
    shutil.copyfile(src=os.path.join(tmpl_dir, 'template.config'), dst=fname)

    # 24 June 2021
    # https://stackoverflow.com/questions/43139174/replacing-words-in-text-file-using-a-dictionary
    for line in fileinput.input(fname, inplace=True):
        line = line.rstrip()
        for key, val in params.items():
            if key in line:
                line = line.replace('@' + key + '@', str(val))
        print(line)


def run_job(config):
    # 24 June 2021
    # https://stackoverflow.com/questions/37058013/how-to-run-a-background-process-and-do-not-wait
    proc = subprocess.Popen(args=['epic', '--config', config],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    print("Job with ID", proc.pid, "submitted.")
    return proc


def wait_running_processes(procs):
    while procs:
        for i, proc in enumerate(procs):
            # 24 June 2021
            # https://stackoverflow.com/questions/43274476/is-there-a-way-to-check-if-a-subprocess-is-still-running
            # https://docs.python.org/3/library/subprocess.html#frequently-used-arguments
            poll = proc.poll()
            if not poll is None:
                if proc.returncode == 0:
                    print("Job with ID", proc.pid, "finished successfully.")
                else:
                    print("Job with ID", proc.pid, "crashed.")
                    exit(-1)
                del procs[i]


try:
    parser = argparse.ArgumentParser(description="Run EPIC simulations.")

    required = parser.add_argument_group('required arguments')

    required.add_argument("--filenames",
                          type=str,
                          nargs='+',
                          required=True,
                          help="list of json files")

    parser.add_argument("--run_dir",
                        type=str,
                        required=False,
                        default=os.getcwd(),
                        help="running directory (default: current working directory)")


    if not '--filenames' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    # 24 June 2021
    # https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    if shutil.which('epic') is None:
        raise RuntimeError('EPIC executable not found in PATH environment variable')

    for fname in args.filenames:
        setup = load_json(fname)
        name = setup['name']


        constants = setup['constants']

        print('Running', name, 'in directory:', args.run_dir)

        for key, values in setup['parameters'].items():

            # the process list
            processes = []

            params = reset_parameters()

            if not key.upper() in params.keys():
                raise KeyError("Parameter '" + key.upper() + "' not available.")

            for kk in constants.keys():
                if kk in params.keys():
                    params[kk] = constants[kk]

            num = 1
            for val in list(values):
                cwd = os.path.join(args.run_dir, name)
                if not os.path.isdir(cwd):
                    os.mkdir(cwd)
                os.chdir(cwd)

                # create config file
                basename = name.lower() + '_' + key.lower() + '_' + str(num)
                params['H5_BASENAME'] = "'" + basename + "'"
                config = basename + '.config'
                create_config_file(params, args.run_dir, config)

                # run simulation
                proc = run_job(config)

                # add to process list
                processes.append(proc)

                os.chdir(args.run_dir)
                num += 1

            wait_running_processes(processes)

            # add to database




except Exception as ex:
    print(ex)
