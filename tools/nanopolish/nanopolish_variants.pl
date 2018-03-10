#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use autodie;

use BioX::Seq::Stream;
use Cwd qw/getcwd abs_path/;
use File::Copy qw/copy/;
use File::Temp qw/tempdir/;
use Getopt::Long qw/:config pass_through/;
use List::Util qw/min/;
use threads;
use threads::shared;

my $fn_genome;
my $threads = 1;
my $fn_outfile;
my $fn_consensus;
my $fn_fast5;
my $fn_reads;
my $fn_index;
my $fn_bam;

# remember full command string (with proper binary)

# parse genome filename and add back to arg stack
GetOptions(
    'genome=s'    => \$fn_genome,
    'threads=i'   => \$threads,
    'outfile=s'   => \$fn_outfile,
    'consensus=s' => \$fn_consensus,
    'fast5=s'     => \$fn_fast5,
    'reads=s'     => \$fn_reads,
    'index=s'     => \$fn_index,
    'bam=s'       => \$fn_bam,
);

my $ret;

my $cwd = abs_path( getcwd() );

$fn_genome    = abs_path( $fn_genome    ) if ( defined $fn_genome    );
$fn_outfile   = abs_path( $fn_outfile   ) if ( defined $fn_outfile   );
$fn_consensus = abs_path( $fn_consensus ) if ( defined $fn_consensus );
$fn_fast5     = abs_path( $fn_fast5     ) if ( defined $fn_fast5     );
$fn_reads     = abs_path( $fn_reads     ) if ( defined $fn_reads     );
$fn_index     = abs_path( $fn_index     ) if ( defined $fn_index     );

# BAM filename is already symbolic link, so abs_path() won't work properly.
# What we actually want are new symbolic links to the same targets
my $fn_bai = $fn_bam;
$fn_bai =~ s/\.bam$/.bai/
    or die "Failed to replace extension of BAM index file";
my $bam_tgt = readlink $fn_bam;
my $bai_tgt = readlink $fn_bai;

my $tmpdir = tempdir( CLEANUP => 1 );

chdir $tmpdir;
mkdir 'tmp';
symlink $bam_tgt, 'input.bam';
symlink $bai_tgt, 'input.bai';

my $fn_link = 'reads';

# divide available threads between actual threads and regions
#
# testing suggests minimal speed-up past 4-8 actual threads per region, so use
# remaining threads for running parallel regions
my $n_threads = min( 4, $threads );
my $n_workers = int($threads/$n_threads);


# extract FAST5 files to path where they are expected
# use system 'tar' to transparently and safely handle absolute paths
my $fast5_dir = 'fast5';
mkdir $fast5_dir;
chdir $fast5_dir;
$ret = system(
    'tar',
    '-xf',
    $fn_fast5
);
die "Failed to extract tarball: $!\n"
    if ($ret);
chdir $tmpdir;

symlink( $fn_reads, $fn_link )
    or die "Failed to create symlink";


# index reads
if (defined $fn_index) {
    $ret = system(
        'tar',
        '-xf',
        $fn_index
    );
    die "Failed to extract tarball: $!\n"
        if ($ret);
}
else {
    $ret = system(
        'nanopolish',
        'index',
        '--directory' => $fast5_dir,
        $fn_link,
    );
    die "Failed nanopolish indexing: $!\n"
        if ($ret);
}

my @cmd = @ARGV;
unshift @cmd, 'nanopolish';
push @cmd, '--genome',  $fn_genome;
push @cmd, '--reads',   $fn_link;
push @cmd, '--threads', $n_threads;
push @cmd, '--bam',     'input.bam';

my @regions :shared;

# build region tags to pass to nanopolish
if (-s $fn_genome) { # gracefully handle empty inputs
    my $parser = BioX::Seq::Stream->new($fn_genome);
    while (my $seq = $parser->next_seq) {
        push @regions, join( ':', $seq->id,
            join( '-', 1, length($seq) ),
        );
    }
}

my @workers;
for (1..$n_workers) {
    push @workers, threads->create(\&run);
}

$_->join() for (@workers);

my @fa_files  = glob "tmp/*.fasta";
my @out_files = glob "tmp/*.vcf";

open my $out_cons, '>', $fn_consensus
    or die "Failed to open output consensus: $!";
for (@fa_files) {
    my $parser =  BioX::Seq::Stream->new($_);
    while (my $seq = $parser->next_seq) {
        $seq->id =~ s/^.+\K:\d+-\d+$//; # strip coordinates from ID
        print {$out_cons} $seq->as_fasta;
    }
}
close $out_cons;

# we may need to do extra processing on VCF output
open my $out_vcf, '>', $fn_outfile
    or die "Failed to open output file: $!";
for my $i (0..$#out_files) {
    my $v = $out_files[$i];
    open my $in, '<', $v;
    while (my $line = <$in>) {
        next if ($line =~ /^\s*#/ && $i > 0);
        print {$out_vcf} $line;
    }
    close $in;
}
close $out_vcf;


sub run {

    LOOP:
    while (1) {

        my $tag;

        {
            lock @regions;
            last LOOP if (! scalar @regions);
            $tag = shift @regions;
        }

        my $fn_out  = "tmp/$tag.vcf";
        my $fn_cons = "tmp/$tag.fasta";

        my @cmd_local = @cmd;
        push @cmd_local, '--window', $tag;
        push @cmd_local, '--outfile',   $fn_out;
        push @cmd_local, '--consensus', $fn_cons;

        my $ret = system @cmd_local;

        my $cmd_string = join ' ', @cmd_local;
        die "Non-zero exit value for command: $cmd_string\n"
            if ($ret);

    }

}
