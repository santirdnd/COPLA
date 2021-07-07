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
use File::Basename qw[basename];

my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: get_ani_identity.pl qry_path ref_fofn\n";
    exit;
}

my $fname_qry = $ARGV[0];
my $ref_fofn = $ARGV[1];
my $work_dir = $fname_qry . '_ani/';
my $fname_tmp = $work_dir . basename($fname_qry) . '.tmp';

my $len_qry;
my $seq = '';
my $line;

# Join multifasta contigs into one contig. TODO: Check other strategies (row of Ns, sort by length, ...)
system("mkdir -p $work_dir") or die;
open my $OUT, '>', "$fname_tmp" or die "Unable to write to $fname_tmp: $!\n";
open my $IN, '<', "$fname_qry" or die "Unable to read from $fname_qry: $!\n";
$line = <$IN>;
print $OUT $line;
while ($line = <$IN>) {
    chomp($line);
    $seq .= $line unless ($line =~ /^>/);
}
print $OUT $seq . "\n";
close $IN;
close $OUT;
$len_qry = length($seq);

#my $threads = '70%';
#my $threads_conf = $work_dir . 'parallel_threads.conf';

#system("echo $threads > $threads_conf");
#system("cat $ref_fofn | parallel -j $threads_conf bin/get_ani_identity_job.pl $fname_tmp $len_qry {}");
system("cat $ref_fofn | parallel bin/get_ani_identity_job.pl $fname_tmp $len_qry {}") or die;

opendir my $DIR, $work_dir or die "Unable to open directory $work_dir: $!\n";
my @files = grep(/\.tsv$/, readdir($DIR));
closedir $DIR;

my $fname_tsv = $fname_qry . '.ani.tsv';
my $fname_dat;
my $data;
open $OUT, '>', "$fname_tsv" or die "Unable to write to $fname_tsv: $!\n";
for my $fname (@files) {
    $fname_dat = $work_dir . $fname;
    open $IN, '<', "$fname_dat" or die "Unable to read from $fname_dat: $!\n";
    $data = <$IN>;
    close $IN;

    $line = $fname_qry . "\t" . substr($fname, 0, length($fname) - 4) . "\t" . $data;
    print $OUT $line;
}
close $OUT;

system("rm -r $work_dir") or die;
