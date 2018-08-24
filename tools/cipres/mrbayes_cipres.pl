#!/usr/bin/perl

use strict;
use warnings;
use Bio::CIPRES;
use Getopt::Long;
use File::Copy qw/copy/;

use constant TOOL_ID => 'MRBAYES_XSEDE';

my $input;
my $cmd;
my $runtime = 0.5;
my $user = 'unknown';
my $name = 'mrbayes';

GetOptions(
    'input=s'   => \$input,
    'cmd=s'     => \$cmd,
    'runtime=f' => \$runtime,
    'user=s'    => \$user,
    'name=s'    => \$name,
) or die "Error parsing command line: $@\n";

die "Must specify input NEXUS file\n"
    if (! defined $input);

my $u = Bio::CIPRES->new(
    conf     => "$ENV{HOME}/.cipres",
    eu       => $user,
    eu_email => $user,
);

my $job_id = sprintf("%08d", rand(100000000)) . "_$user";

my $submit;

open my $in, '<', $input
    or die "Error opening input file: $@\n";
$submit .= $_ while (<$in>);
close $in;

if (defined $cmd) {
    
    die "Command file specified but input file aleady contains mrbayes block\n"
        if ($submit =~ /^\s*begin mrbayes;/im);
    my $cmd_blk = '';
    open my $in, '<', $cmd;
    $cmd_blk .= $_ while (<$in>);
    close $in;
    die "Command file specified but doesn't contain mrbayes block\n"
        if ($cmd_blk !~ /^\s*begin mrbayes;/im);

    # Print command block to STDOUT so user knows exactly what was run
    print "--mrbayes command block------------------------------------\n";
    print $cmd_blk;
    print "-----------------------------------------------------------\n";

    $submit .= $cmd_blk;
}

my ($n_runs, $n_chains) = parse($submit);

my @params = (
    'tool'                      => TOOL_ID,
    'input.infile_'             => $submit,
    'vparam.runtime_'           => $runtime,
    'vparam.mrbayesblockquery_' => 1,
    'vparam.nruns_specified_'   => $n_runs,
    'vparam.nchains_specified_' => $n_chains,
    'metadata.clientJobId'      => $job_id,
    'metadata.clientJobName'    => "Galaxy Job $job_id",
    'metadata.clientToolName'   => TOOL_ID,
    'metadata.statusEmail'      => 'false',
);

print "Submitting job $job_id to CIPRES\n";
my $job = $u->submit_job( @params );

$job->wait() or die "Error waiting for job to finish: $@";

# Pass through STDOUT/STDERR
my $o = $job->stdout;
my $e = $job->stderr;

print "\n--CIPRES STDOUT----------------------------------------\n\n";
print $o;
warn "$e\n";

# Check exit code
my $exit_code = $job->exit_code;
exit($exit_code) if ($exit_code != 0);

mkdir 'output' or die "Error creating output directory: $@ $!\n";
my $out_regex = qr/^infile\.nex.*\.(con\.tre|vstat|lstat|tstat|pstat|parts|trprobs|mcmc)$/;

for my $file ($job->outputs) {
    next if ($file->group ne 'ALL_FILES');
    next if ($file->name !~ /$out_regex/);
    my $new_name = $file->name;
    $new_name =~ s/^infile\.nex/$name/;
    $file->download( out => "output/$new_name" );
}

exit;

sub parse {

    my ($submit) = @_;

    my $c_runs;
    my $c_chains;
    my $n_runs   = 1;
    my $n_chains = 1;

    while ($submit =~ /^\s*mcmcp?\s+.*nruns=(\d+)/img) {
        ++$c_runs;
        $n_runs = $1;
    }
    die "Too many nruns specifications found!\n" if ($c_runs > 1);

    pos($submit) = 0;

    while ($submit =~ /^\s*mcmcp?\s+.*nchains=(\d+)/img) {
        ++$c_chains;
        $n_chains = $1;
    }
    die "Too many nchains specifications found!\n" if ($c_chains > 1);

    return ($n_runs, $n_chains);

}
