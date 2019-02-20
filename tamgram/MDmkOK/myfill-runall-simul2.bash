#!/bin/bash

# M/M/K

# input values

# options="3 m f 1 1 0.5 0.5 -1 out 0 1 1 1 1 1 1 1 0"

K="1 20"  # initial and final number of servers

# alloptions=('0.0004400677458 0.01847404397 0.000009237021984 0.99' '0.00002304646204 0.01847404397 0.000009237021984 0.8' '0.00001153762426 0.01847404397 0.000009237021984 0.6' '0.0007041083933 0.01847404397 0.000009237021984 0.99' '0.00003687433926 0.01847404397 0.000009237021984 0.8' '0.00001846019882 0.01847404397 0.000009237021984 0.6' '0.000836128717 0.01847404397 0.000009237021984 0.99' '0.00004378827787 0.01847404397 0.000009237021984 0.8' '0.0000219214861 0.01847404397 0.000009237021984 0.6')


alloptions=('0.00003687433926 0.01847404397 0.000009237021984 0.8' '0.00001846019882 0.01847404397 0.000009237021984 0.6')

# alloptions=('0.01  0.003  0.0001 0.99' '0.01  0.003  0.0001 0.99')


ncases="0  $((${#alloptions[@]}-1))"


for i in `seq $ncases`
do
    
	      options=${alloptions[i]}


				# main program

				PAC_RATECONST=`echo $options | awk '{print $1}'`
				SERVICE_TIME1CONST=`echo $options | awk '{print $2}'`
				SERVICE_TIME2CONST=`echo $options | awk '{print $3}'`
				ALPHACONST=`echo $options | awk '{print $4}'`



				type tgif >/dev/null 2>&1 || { echo >&2 "I require tgif but it's not installed."; exit 1; }

				  

				optionsnsp=${options//[[:blank:]]/_}
				# mkdir $optionsnsp
				cp template.obj $optionsnsp
				cp Makefile.tgm $optionsnsp/Makefile
				pushd $optionsnsp

				# ./command.bash $options 

				for i in `seq $K`;
				do
					if [ -f out$i.zip ]; then
						echo "file out$i.zip exists"
					else
						echo "file out$i.zip does not exist"
					
					fi
	
				done

				popd

done

