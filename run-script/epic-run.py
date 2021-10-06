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
        "FIELD_FILE": None,
        "FIELD_TOL": None,
        "H5_FIELD_FREQ": None,
        "H5_PARCEL_FREQ": None,
        "H5_FIELD_STATS_FREQ": None,
        "H5_PARCEL_STATS_FREQ": None,
        "H5_BASENAME": None,
        "H5_WRITE_FIELDS": None,
        "H5_WRITE_PARCELS": None,
        "H5_WRITE_FIELD_STATS": None,
        "H5_WRITE_PARCEL_STATS": None,
        "H5_OVERWRITE": None,
        "N_PER_CELL": None,
        "LAMBDA_MAX": None,
        "MIN_VRATIO": None,
        "MAX_VRATIO": None,
        "CORRECTION_ITERS": None,
        "GRADIENT_PREF": None,
        "MAX_COMPRESSION": None,
        "LIMIT": None,
        "ALPHA": None,
        "PRECISE_STOP": None,
    }


def load_json(fname):
    # 24 June 2021
    # https://stackoverflow.com/questions/20199126/reading-json-from-a-file
    with open(fname) as f:
        setup = json.load(f)
    return setup


def dump_json(data, fname):
    # 24 June 2021
    # https://stackoverflow.com/questions/12309269/how-do-i-write-json-data-to-a-file
    with open(fname, "w") as f:
        json.dump(data, f, indent=4)


def create_config_file(params, tconfig, fname):
    # 24 June 2021
    # https://stackoverflow.com/questions/123198/how-can-a-file-be-copied
    shutil.copyfile(src=tconfig, dst=fname)

    # 24 June 2021
    # https://stackoverflow.com/questions/43139174/replacing-words-in-text-file-using-a-dictionary
    for line in fileinput.input(fname, inplace=True):
        line = line.rstrip()
        for key, val in params.items():
            if key in line:
                line = line.replace("@" + key + "@", str(val))
        print(line)


def run_job(config):
    # 24 June 2021
    # https://stackoverflow.com/questions/37058013/how-to-run-a-background-process-and-do-not-wait
    proc = subprocess.Popen(
        args=["epic", "--config", config],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
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

    required = parser.add_argument_group("required arguments")

    required.add_argument(
        "--filenames", type=str, nargs="+", required=True, help="list of json files"
    )

    parser.add_argument(
        "--run_dir",
        type=str,
        required=False,
        default=os.getcwd(),
        help="running directory (absolute path) (default: current working directory)",
    )

    required.add_argument(
        "--tconfig", type=str, required=True, help="template configuration file"
    )

    parser.add_argument(
        "--sequential",
        required=False,
        action="store_true",
        help="run simulations after each other (sequential scan)",
    )

    if (not "--filenames" in sys.argv) or (not "--tconfig" in sys.argv):
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    if not os.path.isfile(args.tconfig):
        raise RuntimeError("Template file not found.")

    if not os.path.isabs(args.run_dir):
        raise RuntimeError("Relative path of running directory.")

    # 24 June 2021
    # https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    if shutil.which("epic") is None:
        raise RuntimeError("EPIC executable not found in PATH environment variable")

    for fname in args.filenames:
        setup = load_json(fname)
        name = setup["name"]

        output = {}

        constants = setup["constants"]

        cwd = os.path.join(args.run_dir, name)
        if not os.path.isdir(cwd):
            os.mkdir(cwd)
        os.chdir(cwd)

        print("Running in directory:", cwd)

        for key, values in setup["parameters"].items():

            key = key.upper()

            # the process list
            processes = []

            params = reset_parameters()

            if not key in params.keys():
                raise KeyError("Parameter '" + key + "' not available.")

            for kk in constants.keys():
                if kk in params.keys():
                    params[kk] = constants[kk]

            num = 1
            for val in list(values):

                # update parameter
                params[key] = val

                # create config file
                basename = name.lower() + "_" + key.lower() + "_" + str(num)
                params["H5_BASENAME"] = "'" + basename + "'"
                config = basename + ".config"
                create_config_file(params, args.tconfig, config)

                # add to database
                output[basename] = params.copy()

                # run simulation
                proc = run_job(config)

                # add to process list
                processes.append(proc)

                if args.sequential:
                    wait_running_processes(processes)

                num += 1

            wait_running_processes(processes)

        dump_json(output, os.path.join("output.json"))

        os.chdir(args.run_dir)

except Exception as ex:
    print(ex)
