#! /bin/bash

DIR=$1

for i in `ls $DIR`; do
    echo "Processing", $DIR/$i
    python scripts/import_ont_model.py -i $DIR/$i -o etc/r9-models/ > new_files.txt
done

for i in `cat new_files.txt`; do
    echo "Dropping", $i
    python scripts/dropmodel.py -i $i
done

rm new_files.txt
