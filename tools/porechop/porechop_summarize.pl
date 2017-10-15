#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

my $flag = 0;
my $table;

while (my $line = <STDIN>) {

    print STDERR $line;

    if ($line =~ /Barcode\s+Reads\s+Bases\s+File/) {
        $flag = 1;
    }
    next if (! $flag);

    $table .= $line;

}

$table =~ s/(^\s+|\e\[\d+m)//gms;
pos($table) = 0;
$table =~ s/[ \t]+/\t/gms;
say $table;

exit;

