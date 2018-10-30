#!perl

my $ori=$/; 
$/='"stLFR_barcode_frequency": {';
my $junk = <>;
$/="}";
chomp($data=<>);
$/=$ori;
chomp($data);
$data =~ s/\n//g;
@a=split(/\t+/,$data);
for($i=1;$i<@a;$i++){
	my @line = split(/[":]/,$a[$i]);
	print "$line[1]\t$line[3]\n";
}
exit;
