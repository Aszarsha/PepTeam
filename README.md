Pepteam
=======
PepTeam is a tool for biomarker discovery through the analysis of screened repertoires of phage displayed libraries that have been sequenced by NGS.  Peptide repertoires are mapped against the proteome of the target organism and functional analysis is performed assuming that small peptides mimic interactions of larger proteins.

# Installation

## System Requirements

PepTeam is a C++ program developped and intended to be used under a GNU/Linux system

* G++ >= 4.8.1
* Boost libraries >= 1.46.1

## Obtaining PepTeam

The source code repository for the version currently in deployment inside the CBIB is available at

* https://github.com/cbib/pepteam

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
### PepteamScoring 

Give a score for each protein depending on peptides mapped
	input :  PepteamProfile output file
	output : 2 files, one in tsv and one in csv showing protein with significant score

	Usage: bin/clt_test.py -i peptides.fastIdx.pepTree.k.mapping.h.profiles

###Pepteam annotation


    	input : tsv file produced by the previous script clt_test.py. this file contains 3 columns, one for protein ID, one for p-value(>0.05) and finally the Z-score
    	output : a csv annotatted file 
    	here is a example of annotated file:


    	"Ensembl proteinId" ,"p-value"     ,"Zscore","Ensembl geneID"    ,EntrezgeneId","associated gene name","description", "wikigene description" 
    	"ENSMUSP00000037233","6.797434e-09","5.6784","ENSMUSG00000039967","30046","Zfp292","zinc finger protein 292 [Source:MGI Symbol;Acc:MGI:1353423]","zinc finger protein 292" 
    	"ENSMUSP00000095766","7.122627e-09","5.6704","ENSMUSG00000039967","30046","Zfp292","zinc finger protein 292 [Source:MGI Symbol;Acc:MGI:1353423]","zinc finger protein 292" 
    
    Usage: pepteamAnnot.pl --input1=database_file --input2=significance.tsv --output=output


    where database_file corresponding to the csv annotated proteome with biomart tool from Ensembl (www.ensembl.org/biomart/martview/3d027c353a68582d5ace7229c16e333e) with this attributes: * Ensembl gene/Protein Id. (Attributes/Genes) * Entrez gene Id. (Attributes/External References) * Associated gene name. (Attributes/Genes) * Description. (Attributes/Genes) * Wikigene description. (Attributes/External References)

# Licence

Copyright (c) 2013,   Hayssam Soueidan (2,3), Thomas Hume (1,2), Benjamin Dartigues (2), and Macha Nikolski (1,2)

(1) LaBRI, Université de Bordeaux, 351 cours de la Libération, 33405 Talence Cedex, France   
(2) CBiB, Université de Bordeaux, 146, rue Léo Saignat 33076 Bordeaux   
(3) INSERM U1035, Université de Bordeaux, 146 rue Léo Saignat, 33076 Bordeaux


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
