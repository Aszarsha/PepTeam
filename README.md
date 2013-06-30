Pepteam
=======

FastIdx
-------
Transform an input multi-fasta file into a .fastIdx index.

    g++ -std=c++1y src/FastIdx.cpp src/PepTree.cpp -o bin/FastIdx src/FastIdx_drv.cpp

PepTree
-------
Transform an input .fastIdx index into a serialized .pepTree.x tree structure representing the set of all windows of size x in the input.

    g++ -std=c++1y src/FastIdx.cpp src/PepTree.cpp -o bin/PepTree src/PepTree_drv.cpp

PepteamMap
----------
Map the first input tree onto the second with a given similarity threshold.

    g++ -std=c++1y src/FastIdx.cpp src/PepTree.cpp -o bin/PepteamMap src/PepteamMap.cpp

General use
-----------

Example : in order to map a peptide repertoire A.txt, each peptide of size 7, onto a MusMusculus.fa proteome database, with a similarity threshold of 0.25, do the following :

	bin/FastIdx -c A.txt
	bin/FastIdx -c MusMusculus.fa
	bin/PepTree -c A.txt.fastIdx 7
	bin/PepTree -c MusMusculus.fa.fastIdx 7
	bin/PepteamMap A.txt.fastIdx.pepTree.7 MusMusculus.fa.fastIdx.pepTree.7 0.25


Happy mapping. :-)
