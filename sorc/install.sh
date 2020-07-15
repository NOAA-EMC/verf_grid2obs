#!/bin/sh
export SRC=$(pwd)
export LOG=${SRC}/install.log
rm ${LOG}
#rm ${SRC}/*.fd/*.o
for i in *.fd/; do 
    echo -ne "Install ${i}\n" | tee -a ${LOG}
    cd $i
    export EXEC=$(basename `grep "CMD.*\=" makefile |awk '{print $NF}'`)
    ls -l $EXEC
    cp -p $EXEC ${SRC}/../exec/
    export err=$?
    if [ "$err" != "0" ]; then
        echo "Failed in ${i}"| tee -a $LOG
        exit 1
    fi
    cd $SRC
    echo -ne "\n\n"| tee -a ${LOG}
done
