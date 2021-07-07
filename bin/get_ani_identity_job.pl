#!/usr/bin/env perl

BEGIN {
    use strict;
    die "Old version of strict module\n" unless strict->VERSION >= 1.0;
    use warnings;
    die "Old version of warnings module\n" unless warnings->VERSION >= 1.1;
    use List::Util;
    die "Old version of List::Util module\n" unless List::Util->VERSION >= 1.41;
    use File::Basename;
    die "Old version of File::Basename module\n" unless File::Basename->VERSION >= 2.85;
}

use strict;
use warnings;
use List::Util qw[min];
use File::Basename qw[basename fileparse];

my $num_args = $#ARGV + 1;
if ($num_args != 3) {
    print "\nUsage: get_ani_identity_job.pl qry_path qry_length ref_path\n";
    exit;
}

my $perc = 0.5;
my $window = 1000;
my $step = 200;

my $fname_qry = $ARGV[0];
my $len_qry = $ARGV[1];
my $fname_ref = $ARGV[2];

my (undef, $work_dir) = fileparse($ARGV[0]);
my $fname_tsv = $work_dir . basename($fname_ref) . '.tsv';

my $len_ref;
my $len_core;
my $threshold;

# Get reference sequence length
my $seq = '';
my $line;
open my $IN, '<', "$fname_ref" or die "Unable to read from $fname_ref: $!\n";
while ($line = <$IN>) {
    chomp($line);
    $seq .= $line unless ($line =~ /^>/);
}
close $IN;
$len_ref = length($seq);

$len_core = $perc * min($len_qry, $len_ref);
if ($len_core <= $window) {
    $threshold = 1;
} else {
    $threshold = 1 + int(($len_core - $window) / $step);
}

system("bin/ani.rb -1 $fname_qry -2 $fname_ref -T $fname_tsv -n $threshold -w $window -s $step -d 3 -t 1 -q > /dev/null 2>&1") or die;
