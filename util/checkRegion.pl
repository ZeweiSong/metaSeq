#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Read rRNA regions from bed file and give each alignment a tag to it.
# Author:            Chao | fangchao@genomics.cn
# ===================================================================
# see detail below

my ($bed,$inf,$out) = @ARGV;

open BED,"< $bed";
open OUT,"> $out";

my (%HS,%C16,%C23);
while(<BED>){
	chomp; my @s = split;
	$HS{$s[0]}{$s[1]}{'end'} = $s[2];
	$HS{$s[0]}{$s[1]}{'srd'} = $s[5];
	if($s[3]=~/16S/){
		$C16{$s[0]}++;	$HS{$s[0]}{$s[1]}{'reg'} = "16S_$C16{$s[0]}";
	}elsif($s[3]=~/23S/){
		$C23{$s[0]}++; $HS{$s[0]}{$s[1]}{'reg'} = "23S_$C23{$s[0]}";
	}else{
		$HS{$s[0]}{$s[1]}{'reg'} = $s[3];
	};
}
close BED;
if($inf=~/m6$/){
	open BLAST, "< $inf";
	while(<BLAST>){
		chomp; my @s = split;
		my ($reg,$match,$frag,$len,$pos,$pct) = ("missing",$s[3],$s[3],"NA",0);
		foreach my $start (sort {$a<=>$b} keys %{$HS{$s[1]}}){
			if($s[8]>=$start-500&&$s[8]<$HS{$s[1]}{$start}{'end'}+500){
				$len = $HS{$s[1]}{$start}{'end'} - $start + 1;
				($match,$frag) = ($s[3],($s[12])?$s[12]:$s[3]);
				$pos = ($HS{$s[1]}{$start}{'srd'} eq "+")?($s[8] - $start + 1):($HS{$s[1]}{$start}{'end'} - $s[8] + 1 - $s[3]);
				$pct = $pos / $len; $st= $start;
				$reg = $HS{$s[1]}{$start}{'reg'}; last;
			}
		}
		my @t = sort {$a<=>$b} ($s[8],$s[9]);
		my $res = sprintf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.4f\n",$s[0],$s[1],$reg,$match,$frag,$t[0],$len,$pos,$pct);
		print OUT $res;
	}
	close BLAST;
}else{
	if($inf=~/bam$/){
		open BAM, "samtools view $inf|";
	}else{
        open BAM, "< $inf";
	}
	while(<BAM>){
		chomp; my @s = split;
		my @f=split("",reverse sprintf("%012b",$s[1]));
		next unless ($f[8] eq 0 && $f[11] eq 0);
		my ($reg,$match,$frag,$len,$pos,$pct) = ("missing",0,0,"NA",0);
		foreach my $start (sort {$a<=>$b} keys %{$HS{$s[2]}}){
			if($s[3]>=$start-100&&$s[3]<$HS{$s[2]}{$start}{'end'}){
				$len = $HS{$s[2]}{$start}{'end'} - $start + 1;
				($match,$frag) = &cigarLen($s[5]);
				$pos = ($HS{$s[2]}{$start}{'srd'} eq "+")?($s[3] - $start + 1):($HS{$s[2]}{$start}{'end'} - $s[3] + 1 - $frag);
				$pct = $pos / $len;
				$reg = $HS{$s[2]}{$start}{'reg'}; last;
			}
		}
		my @d = split("_",$s[0]);
		$d[4] =~ s/multi=//; $d[5] =~ s/len=//;
		my $res = sprintf("%s\t%s\_%s\t%s\t%f\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\n",
		$d[0],$d[1],$d[2],$d[3],$d[4],$d[5],$d[6],$d[8], $s[2],$reg,$match,$frag,$len,$pos,$pct);
		print OUT $res;
	}
	close BAM;
}

close OUT;
exit;

##sub
sub cigarLen{
	my $cigar =shift;
	my ($mLen,$fLen) = (0,0);
	while($cigar){
		$cigar =~ s/^(\d*)([MIDNSHP=X])//;
		my ($mode,$n) = ($2,$1);
		$n ||=1;
		if($mode =~/[MINP=X]/){$mLen += $n;}
		if($mode =~/[HSMINP=X]/){$fLen += $n;}
	}
	return($mLen,$fLen)
}
