#!/usr/bin/env python3

import sys, os
import glob
import subprocess
import shutil
import h5py
import numpy as np

def main():
    tar_file = sys.argv[1]
    out_file = sys.argv[2]
    threads  = sys.argv[3] # currently unused

    extract_fast5(tar_file)

    subprocess.call([
        "nanopolish",
        "extract",
        "--recurse",
        "--fastq",
        "--output", out_file,
        "in_dir" ])

def extract_fast5(fn):

    try:
        in_dir = "in_dir"
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

        # python's tarfile interface does not sanitize file paths within
        # tarballs, which can be a big security risk. GNU tar does sanitize by
        # default, so it's easier/safer here just to call the system tar
        subprocess.call([
            "tar",
            "-xf",
            fn,
            "-C",
            "in_dir"])

        files = glob.glob(
            os.path.join(in_dir, "**", "*.fast5"),
            recursive=True
        )
        if len(files) < 1:
            raise ValueError('No FAST5 files found')

    except OSError as e:
        print("Unexpected error:", e.strerror)
        raise

    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

if __name__ == "__main__" :
    main()
