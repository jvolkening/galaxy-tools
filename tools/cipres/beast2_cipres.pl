#!/usr/bin/perl

use strict;
use warnings;

use Bio::CIPRES;
use XML::LibXML;
use Getopt::Long;
use List::MoreUtils qw/uniq/;

use constant TOOL_ID => 'BEAST2_XSEDE';

my $fn_xml;
my $seed;
my $runtime   = 0.5;
my $user = 'unknown';
my $spec_seed = 0;

GetOptions(
    'in=s'      => \$fn_xml,
    'runtime=f' => \$runtime,
    'spec_seed' => \$spec_seed,
    'seed=i'    => \$seed,
    'user=s'  => \$user,
) or die "Error parsing command line parameters: $@\n";

my ($n_parts, $n_patterns) = parse($fn_xml);

my $u = Bio::CIPRES->new(
    conf     => "$ENV{HOME}/.cipres",
    eu       => $user,
    eu_email => $user,
);

my $has_parts = $n_parts > 1  ? 1 : 0;
my $job_id = sprintf("%08d", rand(100000000)) . "_$user";

my @params = (
    'tool'                    => TOOL_ID,
    'input.infile_'           => [$fn_xml],
    'vparam.runtime_'         => $runtime,
    'vparam.nu_patterns_'     => $n_patterns,
    'vparam.is_partitioned_'  => $has_parts,
    'vparam.spec_seed_'       => $spec_seed,
    'metadata.clientJobId'    => $job_id,
    'metadata.clientJobName'  => "Galaxy Job $job_id",
    'metadata.clientToolName' => TOOL_ID,
    'metadata.statusEmail'    => 'false',
);
push @params, ( 'vparam.nu_partitions_' => $n_parts )
    if ($has_parts);
push @params, ( 'vparam.seed_val_' => $seed )
    if ($spec_seed);

print "Submitting job $job_id to CIPRES\n";
my $job = $u->submit_job( @params );

$job->wait() or die "Error waiting for job to finish: $@";

# Pass through STDOUT/STDERR
my $o = $job->stdout;
my $e = $job->stderr;
print $o;
warn "$e\n";

# Check exit code
my $exit_code = $job->exit_code;
exit($exit_code) if ($exit_code != 0);

for my $file ($job->outputs) {
    next if ($file->group ne 'all_results');
    $file->download( out => $file->name );
}

exit;

sub parse {

    my ($fn_xml) = @_;
    my $dom = XML::LibXML->load_xml( location => $fn_xml )
        or die "Error loading XML: $@\n";
    my %data;
    for my $set ( $dom->findnodes('/beast/data') ) {
        
        my $id = $set->getAttribute('id');
        $data{$id} = $set;

    }

    # Calculate number of variant sites

    my $n_patterns;
    my %data_sets_seen;
    ALN:
    for my $aln ( $dom->getElementsByTagName('alignment') ) {
        my $ref = $aln->getAttribute('idref') || $aln->findvalue('data/@idref') || $aln->getAttribute('data');
        if (defined $ref && length $ref) {
            $ref =~ s/^\@//;
            die "missing dataset for $ref\n" if (! defined $data{$ref});

            # Don't parse duplicated alignments more than once
            next ALN if ( defined $data_sets_seen{ $ref } );
            $data_sets_seen{$ref} = 1;

            $aln = $data{$ref};
        }
            
        my @seqs = map {$_->value} $aln->findnodes('sequence/@value');
        for my $p (0..length($seqs[0])-1) {
            my @bases = uniq map {substr $_, $p, 1} @seqs;
            ++$n_patterns if (scalar(@bases) > 1);
        }
    } 

    # Validate output filenames and determine number of partitions in file
    my $fn_log;
    my $tree_ext;
    my $n_part = 0;
    for my $log ( $dom->findnodes('/beast/run/logger') ) {
        if ($log->getAttribute('id') eq 'tracelog') {
            $fn_log = $log->getAttribute('fileName');
            die "Log filenames must end in '.log' for use by Galaxy\n"
                if ($fn_log !~ /\.log$/);
        }
        elsif ($log->getAttribute('id') =~ /^treelog\.t\:.+$/) {
            my $fn = $log->getAttribute('fileName');
            if ($fn =~ /^.+(\..+)$/) {
                die "Error parsing tree output extension\n"
                    if (defined $tree_ext && $1 ne $tree_ext);
                die "Tree filenames must end in '.trees' for use by Galaxy\n"
                    if ($1 ne '.trees');
            }
            else {
                die "Unexpected mismatch in tree name pattern\n";
            }
            ++$n_part;
        }
    }

    return ($n_part, $n_patterns);

}


