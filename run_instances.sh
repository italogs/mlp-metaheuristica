
make clean
make

NUMEXECS='10'
INSTANCEPATH='instancias'
OUTPUTPATH='outputs'
echo "Instances Path: $INSTANCEPATH"
echo "Output Path: $OUTPUTPATH"
mkdir -p $OUTPUTPATH

echo "" > $OUTPUTPATH/all_outputs.csv

# ./main 'instancias/att48.tsp' 3 

for i in $(ls $INSTANCEPATH); do
	INSTANCEFILE=$INSTANCEPATH/$i;
	OUTPUTFILE=$OUTPUTPATH/$i;
	echo "" > $OUTPUTFILE;
	echo $INSTANCEFILE;
	echo $OUTPUTFILE;
	j=1;
	for (( j = 1; j <=$NUMEXECS; j++ )); do
		echo "<>Execution: $j ----------------------------------</>" >> $OUTPUTFILE;
		echo "sudo ./main $INSTANCEFILE $j >> $OUTPUTFILE;";
		./main $INSTANCEFILE $j >> $OUTPUTFILE;
		echo "<>End Execution: $j ------------------------------</>" >> $OUTPUTFILE;
	done
done
###################################################################