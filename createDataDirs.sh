#!/bin/bash

PARENTDIR='GSE193337_RAW'
mkdir $PARENTDIR

for FILE in *.tsv.gz; do
DIRNAME=`echo $FILE |awk -F_ '{print $1}'`; NEWFILE=`echo $FILE |awk -F_ '{print $3}'`;
    if [ -d "$PARENTDIR/$DIRNAME" ]
    then mv $FILE "$PARENTDIR/$DIRNAME/$NEWFILE"
    else mkdir "$PARENTDIR/$DIRNAME"
        mv $FILE "$PARENTDIR/$DIRNAME/$NEWFILE"
    fi
done

for FILE in *.mtx.gz; do
DIRNAME=`echo $FILE |awk -F_ '{print $1}'`; NEWFILE=`echo $FILE |awk -F_ '{print $3}'`;
    if [ -d "$PARENTDIR/$DIRNAME" ]
    then mv $FILE "$PARENTDIR/$DIRNAME/$NEWFILE"
    fi
done