#!/usr/bin/env python3

import sys, os
import glob
import tarfile
import subprocess
import shutil
import h5py
import numpy as np
from distutils.util import strtobool
from tempfile import mkdtemp

def main():
    tar_file     = sys.argv[1]
    out_file     = sys.argv[2]
    fast5_out    = strtobool( sys.argv[3] )
    demux        = false
    enable_filt  = strtobool( sys.argv[5] )
    hp_correct   = strtobool( sys.argv[6] )
    threads      = sys.argv[7]

    tempdir = mkdtemp()

    (flowcell, kit) = parse_meta(tar_file, tempdir)

    subprocess.call(
        ["guppy_basecaller",
        "--input_path", tempdir,
        "--cpu_threads_per_caller", threads,
        "--save_path", "out_dir",
        "--flowcell", flowcell,
        "--kit", kit,
        "--recursive",
        "--quiet",
		"--disable_pings",
        "--hp_correct", hp_correct,
        "--records_per_fastq", "0" ] +
        ["--qscore_filtering"] * enable_filt +
        ["--fast5_out"] * fast5_out )

    out_path = "out_dir/workspace"
    pass_path = os.path.join( out_path, "pass" )
    if os.path.exists( pass_path ):
        out_path = pass_path
    if demux:
        #check for demuxed albacore output and copy to Galaxy output
        final_dir = "final"
        if not os.path.exists(final_dir):
            os.makedirs(final_dir)
        dirs = glob.glob( os.path.join(out_path, "*") )
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
                if (os.path.isfile(d)):
                    if (d.endswith('.fast5')):
                        bc = os.path.basename( os.path.normpath(d) ) + ".tar.gz"
                        files = [d]
                    else:
                        continue
                else:
                    bc = os.path.basename( os.path.normpath( d ) ) + ".fast5.tar.gz"
                    files = glob.glob( os.path.join( d, "**", "*.fast5"), recursive=True)
                out = os.path.join( final_dir, bc )
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
            files = glob.glob( os.path.join(out_path, "*.fastq") )
            if len(files) != 1:
                raise ValueError('No or multiple FASTQ output files found')
            found_file = files[0]
            shutil.copy(found_file, out_file)
        elif out_fmt == 'fast5':
            #check for single albacore output and copy to Galaxy output
            files = glob.glob( os.path.join(out_path,"**","*.fast5"), recursive=True )
            if len(files) < 1:
                raise ValueError('No FAST5 output files found')
            tar = tarfile.open(out_file, 'w:gz')
            tar.add(out_path)
            tar.close()
        else:
            raise ValueError('Bad output format specified')

    try:
        shutil.rmtree(tempdir)
    except:
        print("Unable to remove temp directory:", sys.exc_info()[0])
        raise


def parse_meta(fn, in_dir):

    try:
        # python's tarfile interface does not sanitize file paths within
        # tarballs, which can be a big security risk. GNU tar does sanitize by
        # default, so it's easier/safer here just to call the system tar
        subprocess.call([
            "tar",
            "--warning=no-unknown-keyword",
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
        test_file = files[0]

        f = h5py.File(test_file,"r")
        #TODO: clean up attribute checking
        try:
            flowcell = f["/UniqueGlobalKey/context_tags"].attrs["flowcell"].upper()
        except:
            try:
                flowcell = f["/UniqueGlobalKey/context_tags"].attrs["flowcell_type"].upper()
            except:
                raise ValueError('No attribute found for flowcell type')
        try:
            kit = f["/UniqueGlobalKey/context_tags"].attrs["sequencing_kit"].upper()
        except:
            raise ValueError('No attribute found for sequencing kit')
            
    except OSError as e:
        print("Unexpected error:", e.strerror)
        raise

    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

    return flowcell, kit

if __name__ == "__main__" :
    main()
