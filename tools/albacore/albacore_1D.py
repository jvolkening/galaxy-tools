#!/usr/bin/env python3

import sys, os
import glob
import tarfile
import subprocess
import shutil
import h5py
import numpy as np
import tarfile
from distutils.util import strtobool

def main():
    tar_file = sys.argv[1]
    out_file = sys.argv[2]
    out_fmt  = sys.argv[3]
    demux    = strtobool( sys.argv[4] )
    threads  = sys.argv[5]

    (flowcell, kit) = parse_meta(tar_file)

    subprocess.call(
        ["read_fast5_basecaller.py",
        "--input", "in_dir",
        "--worker_threads", threads,
        "--save_path", "out_dir",
        "--flowcell", flowcell,
        "--kit", kit,
        "--recursive",
        "--files_per_batch_folder", "0",
        "--output_format", out_fmt,
        "--reads_per_fastq_batch", "999999999" ] +
        ["--barcoding"] * demux )

    if demux:
        #check for demuxed albacore output and copy to Galaxy output
        final_dir = "final"
        if not os.path.exists(final_dir):
            os.makedirs(final_dir)
        dirs = glob.glob("out_dir/workspace/*")
        for d in dirs:

            if out_fmt == 'fastq':
                bc = os.path.basename( os.path.normpath( d ) ) + ".fastq"
                print(d)
                print(bc)
                out = os.path.join( final_dir, bc )
                files = glob.glob( os.path.join( d, "*.fastq") )
                if len(files) != 1:
                    raise ValueError('No or multiple FASTQ output files found')
                found_file = files[0]
                shutil.copy(found_file, out)

            elif out_fmt == 'fast5':
                bc = os.path.basename( os.path.normpath( d ) ) + ".fast5.tar.gz"
                print(d)
                print(bc)
                out = os.path.join( final_dir, bc )
                files = glob.glob( os.path.join( d, "**", "*.fast5"), recursive=True)
                if len(files) < 1:
                    raise ValueError('No FAST5 output files found')
                tar = tarfile.open(out, 'w:gz')
                tar.add( d )
                tar.close()

            else:
                raise ValueError('Bad output format specified')

    else:
        if out_fmt == 'fastq':
            #check for single albacore output and copy to Galaxy output
            files = glob.glob("out_dir/workspace/*.fastq")
            if len(files) != 1:
                raise ValueError('No or multiple FASTQ output files found')
            found_file = files[0]
            shutil.copy(found_file, out_file)
        elif out_fmt == 'fast5':
            #check for single albacore output and copy to Galaxy output
            files = glob.glob("out_dir/workspace/**/*.fast5", recursive=True)
            if len(files) < 1:
                raise ValueError('No FAST5 output files found')
            tar = tarfile.open(out_file, 'w:gz')
            tar.add("out_dir/workspace")
            tar.close()
        else:
            raise ValueError('Bad output format specified')


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
