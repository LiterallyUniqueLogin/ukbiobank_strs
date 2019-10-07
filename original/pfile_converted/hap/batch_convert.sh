#for i in {1..22} ; do
for i in 2 3 9 14 16 20 ; do
	qsub -v "INPUT1=$i" convert.sh &
done

wait
