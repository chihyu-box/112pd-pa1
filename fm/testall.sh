#! /bin/bash

for i in {0..5}
do
	./bin/fm ../input_pa1/input_$i.dat input_$i.out
	./../evaluator/checker ../input_pa1/input_$i.dat input_$i.out
done
