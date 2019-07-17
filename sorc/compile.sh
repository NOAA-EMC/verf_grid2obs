#!/bin/sh
#Floyd Fayton
export SRC=$(pwd)
export LOG=${SRC}/compile.log
rm ${LOG}
rm ${SRC}/*.fd/*.o
for i in *.fd/; do 
    echo -ne "Making ${i}\n" | tee -a ${LOG}
    cd $i
    export EXEC=$(grep "CMD.*\=" makefile |awk '{print $NF}')
    rm $EXEC
    make -f makefile | tee -a $LOG
    export err=$?
    if [ "$err" != "0" ]; then
        echo "Failed in ${i}"| tee -a $LOG
        exit 1
    fi
    mv $EXEC ${SRC}/../exec/
    export err=$?
    if [ "$err" != "0" ]; then
        echo "Failed in ${i}"| tee -a $LOG
        exit 1
    fi
    cd $SRC
    echo -ne "\n\n"| tee -a ${LOG}
done
