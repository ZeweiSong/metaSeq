#!perl

my $ori=$/;
$/='"stLFR_barcode_frequency": {';
my $junk = <>;
$/="}";
chomp($data=<>);
$/=$ori;
chomp($data);
$data =~ s/,\n//g;
@a=split(/\t+/,$data);
my $order = 0;
for($i=1;$i<@a;$i++){
	$order ++;
	my @line = split(/[":]/,$a[$i]);
	printf ("%s\t%d\t%d\n",$line[1],$line[3],$order);
}
exit;
