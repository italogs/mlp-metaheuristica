sudo echo 0 > /sys/devices/system/cpu/cpu4/online
sudo echo 0 > /sys/devices/system/cpu/cpu5/online
sudo echo 0 > /sys/devices/system/cpu/cpu6/online
sudo echo 0 > /sys/devices/system/cpu/cpu7/online

make

NUMEXECS='10'
INSTANCEPATH='instancias2'
OUTPUTPATH='outputs'
echo "Instances Path: $INSTANCEPATH"
echo "Output Path: $OUTPUTPATH"
mkdir -p $OUTPUTPATH


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
		sleep 1;
	done
done
###################################################################