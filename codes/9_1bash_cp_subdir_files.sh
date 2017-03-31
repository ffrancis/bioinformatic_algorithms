#!/bin/bash



### test
# for e in ./output_sampling_21_25_test/*/			### read all directories in the given path
	# do
	# IFS='/' read -r -a array_21_25 <<< "$e"		### break string by '/' and save resulting words as an array
	# y="${array_21_25[2]}"						### select the 3rd item in array
	# cp -R "./output_sampling_21_25_test/$y/." "./output_sampling_1_20_test/$y/"		### copy all subdirectories and files in a given path to another existing directory
# done



### whole data
for e in ./output_sampling_21_25/*/                                         ### read all directories in the given path
	do
	IFS='/' read -r -a array_21_25 <<< "$e"                                 ### break string by '/' and save resulting words as an array
	y="${array_21_25[2]}"                                                   ### select the 3rd item in array
	cp -R "./output_sampling_21_25/$y/." "./output_sampling_1_20/$y/"		### copy all subdirectories and files in a given path to another existing directory
done
	
