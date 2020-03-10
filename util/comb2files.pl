#!/usr/bin/env perl
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       combine 2 files
# Author:            Chao | fangchao@genomics.cn
# ===================================================================
# see detail below

my ($mode,$in1,$in2,$out) = @ARGV;

if ($mode eq "rRNA&clip"){
	open IN1,"< $in1";
	open IN2,"< $in2";
	open OUT,"> $out";
	my %HS;
	while(<IN1>){
		chomp;
		if($_=~/^>/){
				my @s = split />|:/;
				$HS{$s[3]}=$s[1];
				my $reads = <IN1>;
				print OUT "$_\n$reads";
		}
	}
	close IN1;

	while(<IN2>){
		chomp;
		if($_=~/^>(\S+)$/){
			my $id =$1;
			if(defined $HS{$id}){
				<IN2>;
			}else{
				$reads = <IN2>;
				print OUT "$_\n$reads";
			}
		}
	}
	close IN2;
	close OUT;
	exit;

}elsif($mode eq "ucm"){
	open IN1,"< $in1";
	open IN2,"< $in2";
	my (%HS,%RD,$ct);
	while(<IN1>){
		$_=~/^>(\S+)$/;
		chomp($RD{$1}{'id'}=$_);
		chomp($RD{$1}{'seq'}=<IN1>);
	}
	while(<IN2>){
		my @s=split;
		if($s[0] eq "S"){
			$HS{$s[8]}{$s[8]} = "+";
		}elsif($s[0] eq "H"){
			$HS{$s[9]}{$s[8]} = $s[4];
		}elsif($s[0] eq "C"){
			if($s[2]>9){
				$ct ++;
				open OUT,">$out.$ct.fa";
				foreach my $h (keys %{$HS{$s[8]}}){

					if($HS{$s[8]}{$h} eq "-"){
						my $rd = &exchange("rc",$RD{$h}{'seq'});
						print OUT "$RD{$h}{'id'}\n$rd\n";
					}else{
						print OUT "$RD{$h}{'id'}\n$RD{$h}{'seq'}\n";
					}
				}
				close OUT;
			}
		}
	}
}

#########################
sub exchange{
  my %SH;
  my $mod = shift;
  my $str = shift;
  for(my $p=0;$p<length($str);$p++){
    my $pick = substr($str,$p,1);
    my $comp = ($pick eq "A")?"T":(($pick eq "T")?"A":(($pick eq "G")?"C":(($pick eq "C")?"G":$pick)));
    $SH{$p}{'raw'} = $pick;
    $SH{$p}{'com'} = $comp;
  }
  my $res = "";
  if($mod =~ /r/){
    foreach my $p (sort {$b<=>$a} %SH){
      $res .= ($mode =~/c/)?$SH{$p}{'com'}:$SH{$p}{'raw'};
    }
  }else{
    foreach my $p (sort {$a<=>$b} %SH){
      $res .= ($mode =~/c/)?$SH{$p}{'com'}:$SH{$p}{'raw'};
    }
  }
  return($res);
}
