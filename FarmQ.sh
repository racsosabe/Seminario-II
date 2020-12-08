#!/bin/bash

for((i = 24; i <= 50; i++)); do
	echo "Farming Q from $i"
	python3 LLL.py $i > "Q$i.txt"
done
