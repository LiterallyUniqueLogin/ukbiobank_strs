for i in {1..22} ; do
	qsub -v "INPUT1=$i" convert.sh &
done

wait
