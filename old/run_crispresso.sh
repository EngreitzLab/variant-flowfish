source /broad/software/scripts/useuse
reuse -q GridEngine8
reuse -q Python-2.7
reuse -q Java-1.8
source /seq/lincRNA/Ben/VENV_CRISPResso/bin/activate
export PATH=/home/unix/bdoughty/CRISPResso_dependencies/bin:$PATH

echo $1 $2 $3 $4 >> test.txt

CRISPResso -r1 $1 -a $2 -g $3 --cleavage_offset -14 --offset_around_cut_to_plot 10 -o crispresso -n $4 -q 30 -s 15