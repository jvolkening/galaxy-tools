#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use BioX::Seq::Stream;
use List::Util qw/sum/;

my $fn_raw1;
my $fn_raw2;
my $fn_filt1;
my $fn_filt2;
my $fn_bedgraph;
my $fn_qc1;
my $fn_qc2;
my $fn_consensus;
my $fn_out;
my $n_threads = 1;

GetOptions(
    'raw_1=s'     => \$fn_raw1,
    'raw_2=s'     => \$fn_raw2,
    'filt_1=s'    => \$fn_filt1,
    'filt_2=s'    => \$fn_filt2,
    'bedgraph=s'  => \$fn_bedgraph,
    'fastqc_1=s'  => \$fn_qc1,
    'fastqc_2=s'  => \$fn_qc2,
    'consensus=s' => \$fn_consensus,
    'out=s'       => \$fn_out,
    'threads=i'   => \$n_threads,
);


my @counts;
for ($fn_raw1, $fn_raw2, $fn_filt1, $fn_filt2) {
    open my $in, '-|', 'wc', '-l', $_;
    my $ret = <$in>;
    close $in;
    chomp $ret;
    my ($count, $fn) = split ' ', $ret;
    die "line length not multiple of four for $_\n"
        if ($count % 4);
    push @counts, $count/4;	
}


die "raw pair count mismatch\n" if ($counts[0] != $counts[1]);
die "filtered pair count mismatch\n" if ($counts[2] != $counts[3]);

# read fragment stats from STDIN
my @lens;
while (<STDIN>) {
    chomp $_;
    push @lens, $_;
}

my $frag_mean = int( sum(@lens)/scalar(@lens)+0.5 );
my $frag_sd = int( sqrt( sum( map {($_ - $frag_mean)**2} @lens)/(scalar(@lens)-1) )+0.5 );

# extract FastQC data
#warn "extracting FastQC stats...\n";

my @five_nums;
for my $fn ($fn_qc1, $fn_qc2) {
    open my $in, '<', $fn;

    my $in_data = 0;
    my @data;
    LINE:
    while (my $line = <$in>) {
        chomp $line;
        if ($in_data) {
            if ($line =~ /^>>END_MODULE/) {
                last LINE;
            }
            next if ($line =~ /^#/);
            my ($score, $count) = split ' ', $line;
            push @data, [$score,$count];
        }
        elsif ($line =~ /^>>Per sequence quality scores/) {
            $in_data = 1;
        }
    }

    push @five_nums, data_to_5( @data );
}

# Count contigs

my $p = BioX::Seq::Stream->new($fn_consensus);
my %n_contigs;
my @names;
while (my $seq = $p->next_seq) {
   
    my $id = $seq->id;
    push @names, $id;
    while ($seq =~ /[^Nn]+/g) {
        ++$n_contigs{$id};
    }
}	

# Parse assembly depth info
#warn "calculating coverage stats...\n";

open my $in, '<', $fn_bedgraph;

my %cov_5nums;
my %counts;
my $last_end;
my $last_contig;
my $head = <$in>;
while (my $line = <$in>) {
    chomp $line;
    my ($contig,$start,$end,$depth) = split "\t", $line;
    $last_contig //= $contig;
    if ($contig ne $last_contig) {

        my @depths = sort {$a <=> $b} keys %counts;
        my @data;
        for (@depths) {
            push @data, [$_, $counts{$_}];
        }
        $cov_5nums{$last_contig} = data_to_5( @data );
        $last_contig = $contig;
        %counts = ();
        $last_end = undef;
    }

    if (defined($last_end) && $last_end < $start) {
        $counts{0} += $start - $last_end;
    }
    $counts{$depth} += $end - $start;
    $last_end = $end;
}
my @depths = sort {$a <=> $b} keys %counts;
my @data;
for (@depths) {
    push @data, [$_, $counts{$_}];
}
$cov_5nums{$last_contig} = data_to_5( @data );

open my $out, '>', $fn_out;
print {$out} join("\t", 
    '#id',
    'raw_read_pairs',
    'filt_read_pairs',
    'frag_len_mean',
    'frag_len_sd',
    'forward_qual',
    'reverse_qual',
    'n_contigs',
    'coverage_depth',
), "\n";
for my $id (@names) {
    print {$out} join("\t",
        $id,
        $counts[0],
        $counts[2],
        $frag_mean,
        $frag_sd,
        $five_nums[0],
        $five_nums[1],
        $n_contigs{$id},
        $cov_5nums{$id},
    ), "\n";
}
close $out;

exit;


sub data_to_5 {

    my (@data) = @_;
    my $total = sum map {$_->[1]} @data;
    my @five_num;
    my $curr = 0;
    for my $i (0..$#data) {
        $curr += $data[$i]->[1];
        for my $j (0..4) {
            next if (defined $five_num[$j]);
            my $quant = $j/4;
            if ($curr/$total > $quant) {
                $five_num[$j] = $data[$i]->[0];
            }
            elsif ($curr/$total == $quant) {
                $five_num[$j] = $i < $#data
                    ? int(($data[$i]->[0]+$data[$i+1]->[0])/2)
                    : $data[$i]->[0];
            }
        }
    }
    return join('|',@five_num);

}
