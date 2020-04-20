#!/usr/bin/env perl
# Writing sequence into file per bead. Quicker by making index file
#-----------------------------------------------------------------------------
# Author : Chao Fang
# Email  : fangchao@genomics.cn
# Create : Nov 2018
#-----------------------------------------------------------------------------
# see usage below
use strict;
use Tie::File;
use Getopt::Long;

sub usage {
  my $msg = shift;
print <<USAGE;
$msg
usage: perl beadsWrite.pl  --r1 <single fastq> -L <list> -c <cutoff> -s <sufix> -o <output>
       perl beadsWrite.pl  --r1 <read 1 fastq> --r2 <read 2 fastq> -b <barcode> -f <format> -o <output>

default options:
    --r1      single fq file
    -L        barcode frequency list or index
    -c        [int] cutoff
    -s        [1|2] sufix
    -f        format [fq|fa|fq+fa|fq.gz]
    -t        threads to use for compressing
    -k        skip if file exists
    -e        exchange barcode and seq ID (e:put barcode in front of seqID;d:remove barcodes)
    -x        make index only
    -p        prefix. When specificed, output seq into a single file with this prefix.
    -o        output dir
    -v        verbose mode
extract specific barcodes:
    --r2      read2 fq file
    -b        barcode id
    -d        individual barcode id
USAGE
  exit;
};

&usage && exit unless @ARGV;
my ($in1,$in2,$list,$cut,$sfx,$bcode,$idvd,$proc,$skip,$exchange,$index,$pfx,$fmt,$out,$verbose);
GetOptions(
    "r1:s"        => \$in1,
    "r2:s"        => \$in2,
    "L|list:s"    => \$list,
    "c|cutoff:s"  => \$cut,
    "s|sufix:s"   => \$sfx,
    "b|barcode:s" => \$bcode,
    "d|individual:s"=>\$idvd,
    "t|thread=i"  => \$proc,
    "k|skip"      => \$skip,
    "e|exchange:s"  => \$exchange,
    "x|index"     => \$index,
    "p|pfx=s"     => \$pfx,
    "F|format:s"  => \$fmt,
    "o|outdir:s"  => \$out,
    "v|verbose"   => \$verbose,
);
$fmt ||= "fq";
$proc||= 1;
$cut ||= 0;

my (%FILES,%LIST,%REST,%STAT,%IDX,$idxEnable);
my $FHcount=0;
my $firstCatchTime = time();
my $lastCatch = $firstCatchTime;
my @timeS;

if($index){
  &makeindex('make');
}elsif($bcode||$idvd){
  &specifcBarcode($in1,$in2,$bcode,$fmt,$out);
}else{
  unless($list){
    if( -e "$in1.idx"){
      &verbose("index haven't specificed but $in1.idx is found. Use it by default.\n");
    }else{
      &verbose("index haven't found. Generating...\n");
      &makeindex('make');
      &verbose("Use index $in1.idx\n");
    }
    $list = "$in1.idx";
  }
  `mkdir -p $out`;
  &writeAll($in1,$list,$cut,$sfx,$proc,$fmt,$out);
}


exit;

###############################
#
###############################
sub makeindex {
  my $operate = shift;
  if($operate eq 'make'){
    if($in1 =~/.gz$/){
      open IN, "pigz -p $proc -dc $in1|" or die $!;
    }else{
      open IN, "< $in1" or die $!;
    }
    open IDX,"> $in1.idx";
    &verbose("Writing index into $in1.idx\n");
    my $pbc = "";
    while(<IN>){
      my $bc = $1 if $_ =~ /[@\/](\d+_\d+_\d+)\//;
      if ($bc ne $pbc){
        print IDX ($. - 1)."\n" if $. > 1;
        print IDX "$bc\t$.\t";
        $STAT{'b'} ++;
        if($STAT{'b'} % 1000000 == 0 || time() - $lastCatch > 60 ){
          &verbose(sprintf("Generating %8d barcodes ...\n",$STAT{'b'}));
        }
        $pbc = $bc;
      }
      <IN>;<IN>;<IN>;
    }
    print IDX "$.\n";
    close IDX;
    close IN;
  }elsif($operate eq 'add'){
    my ($prv,$bc,$start) = @_;
    print IDX "$prv\n" if $. > 1;
    print IDX "$bc\t$start\t" if $bc;
  }

}

sub specifcBarcode{
  my($in1,$in2,$bcode,$fmt,$out) = @_;
  &verbose("specifcBarcode mode. Start.\n");
  my (%BCS);
  if($bcode =~/\./){
    open BC, "<$bcode" or die "barcode list format error.".$!;
    while(<BC>){
      chomp;
      if($pfx){
        push @{$BCS{'BC'}{0}}, $_;
      }else{
        my @a = split(/\t/,$_);
        push @{$BCS{'BC'}{$a[0]}}, $a[1];
      }
    }
    close BC;
  }elsif($bcode =~/\S+/){
    @{$BCS{0}} = split(",",$bcode);
  }
  if($idvd){
    open DD, "<$idvd" or die "barcode list format error.".$!;
    while(<DD>){
      chomp;
      if($pfx){
        push @{$BCS{'BI'}{0}}, $_;
      }else{
        my @a = split(/\t/,$_);
        push @{$BCS{'BI'}{$a[0]}}, $a[1];
      }
    }
    close DD;
  }
  #my $bnum = @{$BCS{0}};

  if(-e "$in1.idx"){
    $idxEnable = 1;
    open IDX,"< $in1.idx";
    &verbose("Index found. Loading ... ");
    while(<IDX>){
      chomp;
      my @L=split(/\t/,$_);
      $IDX{$L[0]}{'S'} = $L[1];
      $IDX{$L[0]}{'E'} = $L[2];
    }
    close IDX;
    &verbose("Done\n",1);
  }else{
    $idxEnable = 0;
    &verbose("No index found. Pls make one first\n");
    exit;
  }
  my (@R1,@R2);
  tie @R1, "Tie::File", $in1, autochomp => 0 or die $!;
  tie @R2, 'Tie::File', $in2, autochomp => 0 or die $!;
  if(defined $pfx){
    &filesOpen(2,"$out/$pfx");
  }
  my @p;
  foreach my $tag (keys %BCS){
    foreach my $cid (sort {$a<=>$b} keys %{$BCS{$tag}}){
      my $cidName = sprintf("%s%08d",$tag,$cid);
      unless(defined $pfx){
        my $oDir = "$out/$cidName";
        if (-e "$oDir/sort.1.fq" && $skip){
          &verbose("$cidName already exist; skipped\n");
          next;
        }
        system "mkdir -p $oDir" unless -e $oDir;
        &filesOpen(2,"$oDir/sort");
      }
      while(@{$BCS{$tag}{$cid}}){
        # Read one sequence info
        my $bc = shift @{$BCS{$tag}{$cid}};
        my ($S,$E) = ($IDX{$bc}{'S'}-1,$IDX{$bc}{'E'}-1);
        die "[ERR] Can not find $bc in the index. Please check!" if $S < 0;
        $FHcount ++ ;
        &verbose("Find #$FHcount : $cidName <= $bc. Writing\n");

        for(my $i=$S;$i<=$E;$i+=4){
          if($fmt eq "fq"){
            print Q1 &exchange($R1[$i]).$R1[$i+1].$R1[$i+2].$R1[$i+3];
            print Q2 &exchange($R2[$i]).$R2[$i+1].$R2[$i+2].$R2[$i+3];
          }
          if($fmt eq "fa"){
            print A1 $R1[$i].$R1[$i+1];
            print A2 $R2[$i].$R2[$i+1];
          }
          $STAT{'barcodes'} ++;
        }
      }
      if($fmt eq "fq"){ close Q1; close Q2;}
      if($fmt eq "fa"){ close A1; close A2;}
    }
  }
  &verbose("$FHcount beads with  $STAT{'barcodes'} reads in total written. Done.\n")
}

sub writeAll {
  my ($in1,$list,$cut,$sfx,$proc,$fmt,$out) = @_;
  &verbose("WriteAll mode. Start.\nReading barcode list ... ");
  open LST,"<$list" or die $!;
  my $fns=0;
  my $rpb=0;
  while(<LST>){
  	chomp;
  	my @str = split;
    $rpb=($#str==2)?($str[2]-$str[1]+1)/4:$str[1];
  	if($rpb>=$cut){
  		$LIST{$str[0]} = $rpb;
  		$fns++;
  	}else{
      $REST{$str[0]} = $rpb;
  		#last;
  	}
  }
  my $f10 = (int $fns/10);
  $f10 = ($f10==0)?1:$f10;
  close LST;
  &verbose("Done. $fns barcodes to write\n",1);
  # Test whether index file exist

  open IN, "pigz -p $proc -dc $in1|" or die $!;

  if(defined $pfx){
    &filesOpen(0,"$out/$pfx");
  }else{
    &filesOpen(3);
  }

  my(@p);
  my $i=0;
  while(<IN>){
    # Read one sequence info
    chomp;
    my @bc = &getIdBcode($_);
    my $id = &exchange($_);
    my $seq = <IN>;
    my $qua = <IN>.<IN>;
    # Detect barcode

    if (defined $LIST{$bc[1]}){
      if($bc[1] ne $p[1]){
        unless(defined $pfx){
          if($FHcount > 0 ){
            if($fmt =~ /fq/){close QT};
            if($fmt =~ /fa/){close AT};
          }
          `mkdir -p $out/$bc[0]`;
          &filesOpen(1,"$bc[0]/$bc[1]");
        }
        $FHcount ++ ;
        if($FHcount % $f10 == 0){
          @timeS =  &getTime(time,$firstCatchTime);
          print STDERR sprintf("[%02d:%02d:%02d] Find #%04s / %04s beads. writing...\n",
            $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount,$fns);
        }
      }
      if($fmt =~ /fq/){print QT $id.$seq.$qua if fileno(QT)};
      if($fmt =~ /fa/){print AT $id.$seq if fileno(AT)};
      $STAT{$bc[1]} ++;
    }elsif (defined $REST{$bc[1]}){
      if($fmt =~ /fq/){print QX $id.$seq.$qua if fileno(QX);}
      if($fmt =~ /fa/){print AX $id.$seq if fileno(AX);}
    }elsif($bc[1] =~ /0000/){
      if($fmt =~ /fq/){print Q0 $id.$seq.$qua if fileno(Q0);}
      if($fmt =~ /fa/){print A0 $id.$seq if fileno(A0);}
    }else{
      if($fmt =~ /fq/){print Q1 $id.$seq.$qua if fileno(Q1);}
      if($fmt =~ /fa/){print A1 $id.$seq if fileno(A1);}
    }

    if($bc[1] ne $p[1]){
      &makeindex('add',$.-1,$bc[1],$.) if fileno(IDX);
      @p = @bc;
    }
  }
  &makeindex('add',$.) if fileno(IDX);

  close IN;
  if($fmt =~ /fq/){
    close QT;
    close Q0;
    close Q1;
    close QX;
  }
  if($fmt =~ /fa/){
    close AT;
    close A0;
    close A1;
    close AX;
  }

  @timeS =  &getTime(time,$firstCatchTime);
  print STDERR sprintf("[%02d:%02d:%02d] %d beads written. Reporting...\n",
    $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0],$FHcount);

  foreach my $b1 (sort keys %STAT){
    print STDOUT "$b1\t$STAT{$b1}\n";
  }
  #close INF;

  @timeS =  &getTime(time,$firstCatchTime);
  print STDERR sprintf("[%02d:%02d:%02d] All done.\n",
    $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0]);

}


sub getTime {
  my $now = shift;
  my $start = shift;
  my ($diff,$d,$h,$m,$s);
  $diff = $now - $start;
  $h = $diff%(3600*24);
  $d = ($diff -$h)/(3600*24);
  $m = $h%3600;
  $h = ($diff -$m)/3600;
  $s = $m%60;
  $m = ($m - $s)/60;
  return($s,$m,$h,$d);
}

sub verbose{
  return 0 unless $verbose;
  my ($msg,$simple) = @_;
  if($simple){
    print STDERR $msg;
  }else{
    my $thisCatch = time();
    @timeS = &getTime($thisCatch,$firstCatchTime);
    print STDERR sprintf("[%02d:%02d:%02d] $msg", $timeS[3]*24 + $timeS[2],$timeS[1],$timeS[0]);
    $lastCatch = $thisCatch;
  }

}


sub getIdBcode {
  my $id = shift;
  $id =~ /[@\/](\d+)_(\d+)_(\d+)\//;
  my @bcode= ($1,$2,$3);
  my @bc = ($bcode[0],"$bcode[0]_$bcode[1]_$bcode[2]");
  return(@bc);
}

sub exchange {
  my $id = shift;
  #chomp($id);
  if($exchange eq "d"){
    $id =~ s/^([@>]\S+)(\/\d+_\d+_\d+)(\/[12])$/\1\3/;
  }elsif(defined $exchange){
    $id =~ s/^([@>])(\S+)\/(\d+_\d+_\d+)\/(\/[12])$/\1\3\/\2\/\4/;
  }
  return($id);
}

sub filesOpen{
  my $num = shift;
  if($num == 3){
    if($fmt =~ /fq.gz/){
      unless($skip && -e "$out/0000.$sfx.fq.gz"){ open Q0, "|pigz -p $proc > $out/0000.$sfx.fq.gz" or die $! };
      unless($skip && -e "$out/once.$sfx.fq.gz"){ open Q1, "|pigz -p $proc > $out/once.$sfx.fq.gz" or die $! };
      unless($skip && -e "$out/xxxx.$sfx.fq.gz"){ open QX, "|pigz -p $proc > $out/xxxx.$sfx.fq.gz" or die $! };
    }elsif($fmt =~ /fq/){
      unless($skip && -e "$out/0000.$sfx.fq"){ open Q0, "> $out/0000.$sfx.fq" or die $! };
      unless($skip && -e "$out/once.$sfx.fq"){ open Q1, "> $out/once.$sfx.fq" or die $! };
      unless($skip && -e "$out/xxxx.$sfx.fq"){ open QX, "> $out/xxxx.$sfx.fq" or die $! };
    }
    if($fmt =~ /fa.gz/){
      unless($skip && -e "$out/0000.$sfx.fa.gz"){ open A0, "|pigz -p $proc > $out/0000.$sfx.fa.gz" or die $! };
      unless($skip && -e "$out/once.$sfx.fa.gz"){ open A1, "|pigz -p $proc > $out/once.$sfx.fa.gz" or die $! };
      unless($skip && -e "$out/xxxx.$sfx.fa.gz"){ open AX, "|pigz -p $proc > $out/xxxx.$sfx.fa.gz" or die $! };
    }elsif($fmt =~ /fa/){
      unless($skip && -e "$out/0000.$sfx.fa"){ open A0, "> $out/0000.$sfx.fa" or die $! };
      unless($skip && -e "$out/once.$sfx.fa"){ open A1, "> $out/once.$sfx.fa" or die $! };
      unless($skip && -e "$out/xxxx.$sfx.fa"){ open AX, "> $out/xxxx.$sfx.fa" or die $! };
    }
  }elsif($num == 1){
    my $subDir = shift;
    if($fmt =~ /fq.gz/){
      unless($skip && -e "$out/$subDir.$sfx.fq.gz"){ open QT,"|pigz -p $proc > $out/$subDir.$sfx.fq.gz" or die $! };
    }elsif($fmt =~ /fq/){
      unless($skip && -e "$out/$subDir.$sfx.fq"){ open QT,"> $out/$subDir.$sfx.fq" or die $! };
    }
    if($fmt =~ /fa.gz/){
      unless($skip && -e "$out/$subDir.$sfx.fa.gz"){ open AT,"|pigz -p $proc > $out/$subDir.$sfx.fa.gz" or die $! };
    }elsif($fmt =~ /fa/){
      unless($skip && -e "$out/$subDir.$sfx.fa"){ open AT,"> $out/$subDir.$sfx.fa" or die $! };
    }
  }elsif($num==0){
    my $sPFX = shift;
    if($fmt =~ /fq.gz/){
      unless($skip && -e "$sPFX.fq.gz"){ open QT,"|pigz -p $proc > $sPFX.fq.gz" or die $! };
    }elsif($fmt =~ /fq/){
      unless($skip && -e "$sPFX.fq"){ open QT,"> $sPFX.fq" or die $! };
    }
    if($fmt =~ /fa.gz/){
      unless($skip && -e "$sPFX.fa.gz"){ open AT,"|pigz -p $proc > $sPFX.fa.gz" or die $! };
    }elsif($fmt =~ /fa/){
      unless($skip && -e "$sPFX.fa"){ open AT,"> $sPFX.fa" or die $! };
    }
  }elsif($num==2){
    my $sPFX = shift;
    if($fmt =~ /fq.gz/){
      unless($skip && -e "$sPFX.1.fq.gz"){ open Q1,"|pigz -p $proc > $sPFX.1.fq.gz" or die $! };
      unless($skip && -e "$sPFX.2.fq.gz"){ open Q2,"|pigz -p $proc > $sPFX.2.fq.gz" or die $! };
    }elsif($fmt =~ /fq/){
      unless($skip && -e "$sPFX.1.fq"){ open Q1,"> $sPFX.1.fq" or die $! };
      unless($skip && -e "$sPFX.2.fq"){ open Q2,"> $sPFX.2.fq" or die $! };
    }
    if($fmt =~ /fa.gz/){
      unless($skip && -e "$sPFX.1.fa.gz"){ open A1,"|pigz -p $proc > $sPFX.1.fa.gz" or die $! };
      unless($skip && -e "$sPFX.2.fa.gz"){ open A2,"|pigz -p $proc > $sPFX.2.fa.gz" or die $! };
    }elsif($fmt =~ /fa/){
      unless($skip && -e "$sPFX.1.fa"){ open A1,"> $sPFX.1.fa" or die $! };
      unless($skip && -e "$sPFX.2.fa"){ open A2,"> $sPFX.2.fa" or die $! };
    }
  }
}
