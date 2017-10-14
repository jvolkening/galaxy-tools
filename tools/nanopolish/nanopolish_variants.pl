#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use Archive::Tar;
use Cwd qw/getcwd abs_path/;
use File::Copy qw/copy/;
use Getopt::Long qw/:config pass_through/;
use threads;
use threads::shared;
use BioX::Seq::Stream;

my $fn_genome;
my $threads = 1;
my $fn_outfile;
my $fn_consensus;
my $fn_fast5;
my $fn_reads;

# remember full command string (with proper binary)

# parse genome filename and add back to arg stack
GetOptions(
    'genome=s'    => \$fn_genome,
    'threads=i'   => \$threads,
    'outfile=s'   => \$fn_outfile,
    'consensus=s' => \$fn_consensus,
    'fast5=s'     => \$fn_fast5,
    'reads=s'     => \$fn_reads,
);

my $tmp_dir = 'tmp_dir';
mkdir $tmp_dir;

$fn_fast5 = abs_path($fn_fast5);

# extract FAST5 files to path where they are expected
my $fast5_dir = 'fast5';
mkdir $fast5_dir;
my $cwd = abs_path( getcwd() );
chdir $fast5_dir;
my $tar = Archive::Tar->new();
$tar->read($fn_fast5);
$tar->extract();
say "done extracting $fn_fast5 here";
chdir $cwd;


# index reads
my $ret = system(
    'nanopolish',
    'index',
    '--directory' => $fast5_dir,
    $fn_reads,
);
die "Failed nanopolish indexing: $!\n"
    if ($ret);

my @cmd = @ARGV;
unshift @cmd, 'nanopolish';
push @cmd, '--genome', $fn_genome;
push @cmd, '--reads', $fn_reads;

my @regions :shared;

# build region tags to pass to nanopolish
my $parser = BioX::Seq::Stream->new($fn_genome);
while (my $seq = $parser->next_seq) {
    push @regions, join( ':', $seq->id,
        join( '-', 1, length($seq) ),
    );
}

my @workers;
for (1..$threads) {
    push @workers, threads->create(\&run);
}

$_->join() for (@workers);

my @fa_files  = glob "$tmp_dir/*.fasta";
my @out_files = glob "$tmp_dir/*.vcf";

open my $out_cons, '>', $fn_consensus
    or die "Failed to open output consensus: $!";
for (@fa_files) {
    open my $in, '<', $_;
    while (my $line = <$in>) {
        print {$out_cons} $line;
    }
    close $in;
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

        my $fn_out  = "$tmp_dir/$tag.vcf";
        my $fn_cons = "$tmp_dir/$tag.fasta";

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
