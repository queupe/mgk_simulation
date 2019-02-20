all:
	gcc createmodel1.c -Wall -g -lm -op2pmodelgenerator
	mv p2pmodelgenerator ~/bin
	g++ ergodicity.cpp -Wall -g -lm -oergodicity
	mv ergodicity ~/bin
clean:
	mv template.obj baktemplate.bak
	rm template* -rf
	mv baktemplate.bak template.obj
	rm out* -rf
	rm file* -rf
	rm *~
	rm err *bak
