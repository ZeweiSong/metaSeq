# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Format fastq file to fit athena programe
# Author:            Chao | fangchao@genomics.cn
# Version:           V0.1
# Last modified:    19 Dec 2018 (since 01 Dec 2018)
# ===================================================================
# see detail below
use strict;
use Getopt::Long;

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage:
  $0 -i input -o output
    -1  fastq r1 file
    -2  fastq r2 file
    -o  output filename (Do not compress)
    -s  sample name (default is '1')
    -k  [T|F] Whether to skip barcodes-undetected reads, default is T.
    -v  verbose
    -h  show help info
USAGE
}
&usage && exit unless @ARGV;

my ($inf1,$inf2,$out,$sampleName,$skip0,$verbose,$help);
GetOptions(
  "1=s" => \$inf1,
  "2:s" => \$inf2,
  "s:s" => \$sampleName,
  "k:s" => \$skip0,
  "o=s" => \$out,
  "v" => \$verbose,
  "h|help|?" => \$help,
);
&usage && exit if $help;

$sampleName ||= 1;
$skip0 ||= "T";
if(not defined $inf1){
  open IN1, "<-" or die $!;
}elsif($inf1 =~ /.gz$/){
  open IN1, "pigz -dc $inf1|" or die $!;
}else{
  open IN1, "<$inf1" or die $!;
}
if($inf2){
  if($inf2 =~ /.gz$/){
    open IN2, "pigz -dc $inf2|" or die $!;
  }else{
    open IN2, "<$inf2" or die $!;
  }
}
open OUT, ($out)?">$out":">STDOUT" or die $!;

# Main
&verbose("[log] Start ... ");
while(<IN1>){
  my ($chunk1,$chunk2);
  my @bc = &getIdBarcode($_);
  my $newID = "$bc[0]/$bc[2]\tBC:Z:$bc[1]-$sampleName\n";
  $chunk1 = $newID.<IN1>.<IN1>.<IN1>;
  if($inf2){
    <IN2>; $chunk2 = $newID.<IN2>.<IN2>.<IN2>;
  }
  next if ($skip0 eq "T" && $bc[1] =~ /0000/);
  print OUT $chunk1.$chunk2;
}

close IN1;
close IN2 if $inf2;
close OUT;
# Main end

&verbose("done!\n");

exit;

sub verbose{
  my $msg = shift;
  print STDERR $msg if $verbose;
}

sub getIdBarcode {
  my $id = shift;
  $id =~ /(^@\S+)\/(\d+)_(\d+)_(\d+)\/([12])/;
  my @bcode= ($2,$3,$4);
  my @bc = ($1,"$bcode[0]_$bcode[1]_$bcode[2]",$5);
  return(@bc);
}
