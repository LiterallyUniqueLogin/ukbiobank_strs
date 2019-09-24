#for i in {1..22} ; do
for i in 5 22 ; do
	qsub -v "INPUT1=$i" convertToPfile.sh &
done

wait
