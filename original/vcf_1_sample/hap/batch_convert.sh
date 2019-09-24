#for i in {1..22} ; do
for i in 12 14 15 4 ; do
	qsub -v "INPUT1=$i" convert.sh &
done

wait
