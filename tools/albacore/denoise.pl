#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use File::Basename qw/basename/;
use File::Copy qw/copy/;
use Getopt::Long;
use List::Util qw/sum/;

my $fn_table;
my @reads;
my @names;
my $min_score = 80;
my $min_frac  = 0.05;
my $rm_unclass = 0;
my $n_keep;
my $fn_summary;

use constant BARCODE => 14;
use constant SCORE   => 15;

my $unclass_tag = 'unclassified';

GetOptions(
    'input=s' => \@reads,
    'name=s'  => \@names,
    'table=s' => \$fn_table,
    'min_score=f' => \$min_score,
    'min_frac=f'  => \$min_frac,
    'remove_unclassified' => \$rm_unclass,
    'n_keep=i'    => \$n_keep,
    'summary=s'   => \$fn_summary,
);

die "Table not found\n"
    if (! -r $fn_table);

open my $tsv, '<', $fn_table;
my $head = <$tsv>;
chomp $head;
my @fields = split "\t", $head;

die "unexpected field order"
    if ($fields[BARCODE] ne 'barcode_arrangement');
die "unexpected field order"
    if ($fields[SCORE] ne 'barcode_score');

my %counts;
my %sums;

while (my $line = <$tsv>) {
    chomp $line;
    my @fields = split "\t", $line;
    my $bin   = $fields[BARCODE];
    my $score = $fields[SCORE];
    $counts{$bin} += 1;
    $sums{$bin}   += $score;
}

if ($rm_unclass) {
    delete $counts{$unclass_tag};
    delete $sums{$unclass_tag};
}

my @keys = sort {$counts{$b} <=> $counts{$a}} keys %counts;

my %scores;

my %status;

@status{ @keys } = ('discarded') x scalar(@keys);

for (@keys) {
    $scores{$_} = $sums{$_}/$counts{$_};
}
my $sum_count  = sum(values %counts);
my $mean_count = $sum_count / scalar(values %counts);

if (defined $n_keep) {
    
    my %rank_scores;
    @rank_scores{ @keys } = map {$scores{$_} * $counts{$_}/$sum_count} @keys;

    @keys = sort {$rank_scores{$b} <=> $rank_scores{$a}} keys %rank_scores;

    @keys = @keys[0..$n_keep-1];

}

else {

    @keys = grep {$scores{$_} >= $min_score} @keys;
    @keys = grep {$counts{$_} >= $mean_count*$min_frac} @keys;

}

@status{ @keys } = ('kept') x scalar(@keys);

@keys = sort {$counts{$b} <=> $counts{$a}} keys %status;

# print summary table
open my $summary, '>', $fn_summary
    or die "Failed to open summary: $!\n";
say {$summary} join "\t",
    '#bin',
    'n_reads',
    'avg_score',
    'status',
;
for (@keys) {
    say {$summary} join "\t",
        $_,
        $counts{$_},
        sprintf("%0.3f",$scores{$_}),
        $status{$_},
    ;
}
close $summary;

mkdir('outputs');

for my $i (0..$#reads) {

    my $fn = $reads[$i];
    my $name = $names[$i];

    my $base = basename($name);
    my $bin = $base;
    $bin =~ s/\.(fastq|fq|fast5\.tar\.gz)$//i;
    if (defined $status{$bin} && $status{$bin} eq 'kept') {
        copy( $fn, "outputs/$base" );
    }

}
