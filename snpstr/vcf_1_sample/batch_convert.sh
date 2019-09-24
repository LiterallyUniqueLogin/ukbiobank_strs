#for i in {1..22} ; do
for i in 16 ; do
	qsub -v "INPUT1=$i" convert.sh &
done

wait
