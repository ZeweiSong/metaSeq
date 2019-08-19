# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Purpose:           CIRCOS visualization
# Parameters:        1
#                    ---------------------------
#                    bead cluster dir
#                    ---------------------------
# Author:            fangchao@genomics.cn
# Version:           V0.11
# Last modified:     19 Feb 2019 (since 19 Jan 2019)
# ===================================================================
#
if [ ! $1 ];
then
  echo -e "Usage:\n  circos.sh <zoom> <cutoff> <alignDir> <refID> <output Dir> [force]" && exit
fi

zoom=$1
cutoff=$2
aDIR=$3
refID=$4
oDIR=$5
force=$6

#aDIR=$BCDIR/align
bDIR=`dirname $0`
conf="$bDIR/template/circos.conf"

suffix="i$cutoff.x$zoom"

echo "####################################"
echo [CC] info: Alignment dir: $aDIR
echo [CC] info: Output CC dir: $oDIR/circos.$suffix
echo [CC] info: link pos zoom: $zoom
echo [CC] info: circos config: $conf
echo [CC] info: refer ID info: $refID
echo "#### CIRCOS visualization start ####"


echo
if [[ -d "$oDIR/data" && -f "$oDIR/circos.conf" && -z $force ]];
then
  echo "[BC] Direcotry $oDIR already exist."
else
  echo "[BC] Making dir and init."
  mkdir -p $oDIR
  mkdir -p $oDIR/data
  cp ${conf/circos/ticks} $oDIR/
  cp ${conf/circos/ideogram} $oDIR/
  cp ${conf/circos/housekeeping} $oDIR/
  sed 's/DEFAULTSUF/'$suffix'/g' $conf > $oDIR/circos.conf
fi

echo
karyoREF="$oDIR/data/karyotype.REF.$suffix.tsv"
linkS2R="$oDIR/data/link.scaf2REF.$suffix.tsv"
if [[ -f $karyoREF && -z $force ]];
then
  echo "[BC] karyotype.REF.tsv exists. Skiped (add \$3 to force re-run)"
else
  echo "[BC] Getting karyotype.REF.tsv."
  metabbq findRegion -r $refID -i $aDIR/scaf2ref.blast6 -o $karyoREF \
  -L $linkS2R -c $cutoff -x $zoom
fi

echo
karyoScaf="$oDIR/data/karyotype.scaf.tsv"
if [[ -f $karyoScaf && -z $force ]];
then
  echo "[BC] karyotype.scaf.tsv exists. Skiped (add \$3 to force re-run)"
else
  echo "[BC] Getting karyotype.REF.tsv."
  grep "^@SQ" $aDIR/read2scaf.sam| \
  awk -v z=$zoom '{sub("SN:","",$2);split($2,arr,"_");sub("LN:","",$3); \
  print "chr - "$2"\t"arr[2]"\t0\t"z*($3+50)"\tgrey5"}' >  $karyoScaf
fi

echo
linkS2S="$oDIR/data/link.scaf2scaf.$suffix.tsv"
if [[ -f $linkS2S && -z $force ]];
then
  echo "[BC] link.scaf2scaf.$suffix.tsv exists. Skiped (add \$3 to force re-run)"
else
  echo "[BC] Getting link.scaf2scaf.$suffix.tsv"
  perl -e 'open IN,"<$ARGV[0]";$z=@ARGV[1];while(<IN>){next if $_ =~ /^@/; @a=split;
  next if $a[6] eq "="; next if $a[5] ne "100M"; next if $a[2] eq $prv or $a[2] eq "*";
  $s=$z*$a[3];$e=$z*$a[7];
  print "$a[2]\t$s\t".($s+$z*100)."\t$a[6]\t$e\t".($e+$z*100)."\n";
  $prv=$a[6]}' $aDIR/read2scaf.sam $zoom |sort|uniq > $linkS2S
fi

echo
stack2S="$oDIR/data/stack.read2scaf.$suffix.tsv"
touch $stack2S # skip this step by now

if [[ -f $stack2S && -z $force ]];
then
  echo "[BC] stack.read2scaf.$suffix.tsv. Skiped (add \$3 to force re-run)"
else
  echo "[BC] Getting stack.read2scaf.$suffix.tsv."
  metabbq findRegion -s $aDIR/read2scaf.sam -x $zoom -o $stack2S -r $karyoREF
fi

echo
circosPNG="$oDIR/circos.$suffix.png"
if [[ -f $circosPNG && -z $force ]];
then
  echo "[BC] $circosPNG exists. Skiped (add \$3 to force re-run)"
else
  echo "[BC] Perparing circos configures ..."
#  thick=`grep -E -o "[0-9]+_[0-9]+_[0-9]+" $stack2S|sort|uniq|awk 'END{print int(500/FNR)}'`
#  if [ -z $thick || $thick -gt 50 ]; then thick=50; fi
  thick=5
  cp ${conf/circos/ticks} $oDIR/
  cp ${conf/circos/ideogram} $oDIR/
  sed 's/DEFAULTSUF/'$suffix'/g;s/DEFAULTTHICKNESS/'$thick'/' $conf > $oDIR/circos.conf
  echo "[BC] running circos..."
  cmd="circos -conf $oDIR/circos.conf -outputdir $oDIR -outputfile circos.$suffix"
  echo $cmd && $cmd
fi

echo "Done!"
