Pepteam
=======

# Installation

## System Requirements

PepTeam is a C++ program developped and intended to be used under a GNU/Linux system

* G++ >= 4.8.1
* Boost libraries >= 1.46.1

## Obtaining PepTeam

The source code repository for the version **currently in deployment** inside the CBIB should be available at

* https://github.com/cbib/pepteam

The source code repository for the **current developmental** version is available at

* https://github.com/Aszarsha/PepTeam

## Installing PepTeam

Grab the version that you are interested in using and use Make. For example, for the version currently in deployment do :

	git clone git://github.com/cbib/pepteam
	cd PepTeam
	make

# Usage

## Pipeline

PepTeam is implemented as a multi-steps pipeline

* first preproccess the inputs FASTA files into fastIdx indices files
* create two pepTree acceleration structures from those two fastIdx preprocessed files
* map the two pepTree files one against the other
* extract profile information from the mappings

Example : in order to map a peptide repertoire A.txt, each peptide of size 7, onto a MusMusculus.fa proteome database, with a similarity threshold of 0.25, do the following :

	bin/FastIdx -c A.txt
	bin/FastIdx -c MusMusculus.fa
	bin/PepTree -c A.txt.fastIdx 7
	bin/PepTree -c MusMusculus.fa.fastIdx 7
	bin/PepteamMap A.txt.fastIdx.pepTree.7 MusMusculus.fa.fastIdx.pepTree.7 0.25
	bin/PepteamProfile A.txt.fastIdx.pepTree.7.mapping.0_25 A.txt.fastIdx A.txt.fastIdx.pepTree.7 MusMusculus.fa.fastIdx MusMusculus.fa.fastIdx.pepTree.7

## Details

The detailed usage of each program follows

### FastIdx

Transform an input multi-fasta file into a .fastIdx index.

	Usage: bin/FastIdx -* input-file
	   where * is one of:
	     c -> create the protein index from input fasta file
	     s -> print the number of proteins in the index
	     p -> print the index in human 'interpretable' formatPepTree

### PepTree

Transform an input .fastIdx index into a serialized .pepTree.x tree structure representing the set of all windows of size x in the input.

	Usage: bin/PepTree -c input-FastIdx-file fragments-size
	          create the PepTree file from the input FastIdx file
	   or: bin/PepTree -% pepTree-file
	   where % is one of:
	     d -> print the tree depth
	     v -> print the tree vector in human 'interpretable' format
	     n -> print the tree nodes in human 'interpretable' format
	     l -> print the tree leaves in human 'interpretable' format
	     p -> print the tree leaf positions in human 'interpretable' format

### PepteamMap

Map the first input tree onto the second with a given similarity threshold.

	Usage: bin/PepteamMap pepTree-query-file pepTree-subject-file cutoff-homology

### PepteamProfile

Construct the profiles for each protein of the proteome database with valid mappings

	Usage: bin/PepteamProfile mapping-file query-fastIdx-file query-pepTree-file subject-fastIdx-file subject-pepTree-file

# Licence

Copyright (c) 2013, Thomas Hume (1)

(1) LaBRI, Universite Bordeaux 1, 351 cours de la Liberation, 33405 Talence Cedex, France

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
