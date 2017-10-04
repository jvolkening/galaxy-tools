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
    threads  = sys.argv[3]

    extract_fast5(tar_file)
    with open(out_file, "w") as outfile:
        subprocess.call(["scrappie",
            "raw",
            "--threads", threads,
            "--outformat", "fasta",
            "in_dir" ],
            stdout=outfile )

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
            in_dir])

        files = glob.glob(
            os.path.join(in_dir, "**", "*.fast5"),
            recursive=True
        )
        if len(files) < 1:
            raise ValueError('No FAST5 files found')
        for f in files:
            shutil.copy(f, in_dir)

    except OSError as e:
        print("Unexpected error:", e.strerror)
        raise

    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

if __name__ == "__main__" :
    main()
