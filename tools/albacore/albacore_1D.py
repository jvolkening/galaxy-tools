#!/usr/bin/env python3

import sys, os
import glob
import tarfile
import subprocess
import shutil
import h5py
import numpy as np

def main():
    tar_file = sys.argv[1]
    out_file = sys.argv[2]
    threads  = sys.argv[3]

    (flowcell, kit) = parse_meta(tar_file)

    subprocess.call(["read_fast5_basecaller.py",
        "--input", "in_dir",
        "--worker_threads", threads,
        "--save_path", "out_dir",
        "--flowcell", flowcell,
        "--kit", kit,
        "--recursive",
        "--files_per_batch_folder", "0",
        "--output_format", "fastq",
        "--reads_per_fastq_batch", "999999999" ])

    #check for single albacore output and copy to Galaxy output
    files = glob.glob("out_dir/workspace/*.fastq")
    if len(files) != 1:
        raise ValueError('No or multiple FASTQ output files found')
    found_file = files[0]
    shutil.copy(found_file, out_file)

def parse_meta(fn):

    try:
        in_dir = "in_dir"
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

        tar = tarfile.open(fn, mode='r')
        tar.extractall(path=in_dir)

        files = glob.glob(
            os.path.join(in_dir, "**", "*.fast5"),
            recursive=True
        )
        if len(files) < 1:
            raise ValueError('No FAST5 files found')
        test_file = files[0]

        f = h5py.File(test_file,"r")
        flowcell = f["/UniqueGlobalKey/context_tags"].attrs["flowcell"].upper()
        kit = f["/UniqueGlobalKey/context_tags"].attrs["sequencing_kit"].upper()
    except OSError as e:
        print("Unexpected error:", e.strerror)
        raise

    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

    return flowcell, kit

if __name__ == "__main__" :
    main()
