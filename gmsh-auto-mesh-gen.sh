#!/bin/bash
cwd=$(pwd);
alStart=1
alEnd=2
alSteps=0.02
alDiff=$(echo "scale=3;$alEnd-$alStart" | bc)
NalSteps=$( echo "scale=0;$alDiff / $alSteps" | bc )

sStart=0.75
sEnd=0.1
sSteps=0.05
sDiff=$(echo "scale=3;$sStart-$sEnd" | bc)
NsSteps=$( echo "scale=0;$sDiff / $sSteps" | bc )

for i in $(seq 0 $NalSteps);
do     
	cd $cwd
	alphaValue=$(echo "scale=3;$alStart+$i"*"$alSteps" | bc)
	dirName1="alpha-"$alphaValue
       	if ! [ -d "$dirName1" ]; then
		mkdir $dirName1
	fi
	cd $dirName1
	cwd1=$(pwd)
	for j in $(seq 0 $NsSteps);
		do	
			cd $cwd1
			sValue=$(echo "scale=3;$sStart-$j"*"$sSteps" | bc)
			dirName2="S-0"$sValue
			if ! [ -d "$dirName2" ]; then
				mkdir $dirName2
			fi
			cd $dirName2
			echo $(pwd)
			cp ../../geom.scr .

			old="S=0.75"
			new="S="$sValue
			sed -i "s/$old/$new/" geom.scr

		
			old="a=1"
			new="a="$alphaValue
			sed -i "s/$old/$new/" geom.scr

			gmsh -2 -order 9 geom.scr

		done
done




