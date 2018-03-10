#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use autodie;

use Cwd qw/getcwd abs_path/;
use File::Copy qw/copy/;
use File::Temp qw/tempdir/;
use Getopt::Long qw/:config pass_through/;
use threads;
use threads::shared;
use BioX::Seq::Stream;

my $fn_link = 'reads';

my $fn_fast5;
my $fn_reads;
my $fn_outfile;

# parse genome filename and add back to arg stack
GetOptions(
    'fast5=s'   => \$fn_fast5,
    'reads=s'   => \$fn_reads,
    'outfile=s' => \$fn_outfile,
);

my $ret;

my $cwd = abs_path( getcwd() );

$fn_fast5   = abs_path($fn_fast5);
$fn_reads   = abs_path($fn_reads);
$fn_outfile = abs_path($fn_outfile);

my $tmpdir = tempdir( CLEANUP => 1 );

chdir $tmpdir;

# extract FAST5 files to path where they are expected
my $fast5_dir = 'fast5';
if (-e $fast5_dir) {
    warn "$fast5_dir exists, won't overwrite";
    exit;
}
mkdir $fast5_dir;
chdir $fast5_dir;

# use system 'tar' to transparently and safely handle absolute paths
$ret = system(
    'tar',
    '-xf',
    $fn_fast5
);
die "Failed to extract tarball: $!\n"
    if ($ret);

chdir $tmpdir;

symlink( $fn_reads, $fn_link )
    or die "Failed to create symlink: $@";

# index reads
$ret = system(
    'nanopolish',
    'index',
    '--directory' => $fast5_dir,
    $fn_link,
);
die "Failed nanopolish indexing: $!\n"
    if ($ret);

my @idx_fns = glob "$fn_link.*";
$ret = system(
    'tar',
    '-cf' => $fn_outfile,
    @idx_fns,
);
die "Failed tarball creation: $!\n"
    if ($ret);
