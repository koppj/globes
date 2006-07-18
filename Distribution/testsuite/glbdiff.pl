#!/usr/bin/perl
# PERL script to results from the GLoBES testsuite

$THRESHOLD = 1e-2;

if ($ARGV[0] eq "" || $ARGV[1] eq "")
{
  print "Invalid arguments.\n";
  exit;
}

open(FILE1, $ARGV[0]);
open(FILE2, $ARGV[1]);

$N_MISMATCHES = 0;
$N_LINES      = 0;
$N_UNITARITY  = 0;
while (<FILE1>)
{
  $N_LINES++;
  $LINE1 = $_;
  $LINE2 = <FILE2>;

  while ($LINE1 =~ /Unitarity mismatch/)
  {
    $N_UNITARITY++;
    print "$ARGV[0]: $LINE1";
    $LINE1 = <FILE1>
  }
  while ($LINE2 =~ /Unitarity mismatch/)
  {
    $N_UNITARITY++;
    print "$ARGV[1]: $LINE2";
    $LINE2 = <FILE2>
  }
  
  $LINE1 =~ /^(<.*>)(.*)$/;
  $META1 = $1;
  $DATA1 = $2;
  $LINE2 =~ /^(<.*>)(.*)$/;
  $META2 = $1;
  $DATA2 = $2;

  if ($META1 != $META2)
  {
    print "Meta data mismatch: $META1 <--> $META2\n";
  }

  if ($DATA1 > 0 && abs($DATA2 - $DATA1)/$DATA1 > $THRESHOLD)
  {
    $N_MISMATCHES++;
    print "Data mismatch: $DATA1 <--> $DATA2\n";
    print "Metadata: $META1\n";
    print "\n"
  }
}

print "$N_MISMATCHES out of $N_LINES lines (".100*$N_MISMATCHES/$N_LINES."%) "
      . "contained a relative discrepancy of more than $THRESHOLD.\n";
print "Number of unitarity mismatches: $N_UNITARITY\n";

close(FILE1);
close(FILE2);


