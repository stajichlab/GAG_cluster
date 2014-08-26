#PBS -l nodes=1:ppn=4,walltime=2:00:00 -q js -N fasta -j oe
module load fasta

N=$PBS_ARRAYID
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
  echo "need a command line number or a PBS_ARRAYID"
  exit
fi
CPU=$PBS_NP
if [ ! $CPU ]; then
 CPU=1
fi
INFILE=GAGcluster.fas
OUT=results
DB=proteomes
LIST=list
file=`head -n $N $LIST | tail -n 1`
echo $file
prefix=`basename $file .aa.fasta`
outfile=GAGcluster-vs-$prefix.FASTA.tab
if [ ! -f $OUT/$outfile ]; then
 fasta36 -m 8 -d 0 -b 5 -T $CPU -E 1e-5 $INFILE $file > $OUT/$outfile
 echo "fasta36 -m 8 -d 0 -b 5 -T $CPU -E 1e-5 $INFILE $file > $OUT/$outfile"
fi
