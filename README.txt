PLEK2: a novel method for predicting lncRNA and mRNA based on sequence intrinsic features and Coding-Net model (Upgraded version of PLEK)
===========================================================
PLEK: predictor of long non-coding RNAs and mRNAs based on k-mer scheme (PLEK download | SourceForge.net)
===========================================================


INSTALLATION
-------------
LncRNA participates in many important regulatory activities of organisms. Its biological structure is similar to messenger RNA (mRNA), in order to distinguish between lncRNA and mRNA (messenger RNA) transcripts more quickly and accurately, we upgraded the alignment-free PLEK to PLEK2.


Requirements
------------
+ [Linux]
+ [Python version > = 3.8.5] (http://www.python.org/)
+ [numpy version >= 1.2.2] (http://www.numpy.org/)
$ pip install numpy
+ [pandas version >= 1.3.3] (http://pandas.pydata.org/)
$ pip install pandas
+ [keras version = 2.4.3] 
$ pip install keras==2.4.3
+ [bio version >= 1.3.2]
$ pip install bio
+ [os version]
$ pip install os
+ [itertools]
$ pip install itertools
+ [regex]
$ pip install regex
--------------------------------------------------



Steps:
1.Download PLEK2_model_v3.tar.gz from https://sourceforge.net/projects/plek2/files/ and decompress it.
$ tar zvxf PLEK2_model_v3
2. Compile PLEK2_model_v3
$ cd PLEK2_model_v3
3. decompress Coding_Net_kmer6_orf.h5.bz2 model
$ bunzip2 Coding_Net_kmer6_orf.h5.bz2 
4.  decompress Coding_Net_kmer6_orf_Arabidopsis.h5.bz2 model
$ bunzip2 Coding_Net_kmer6_orf_Arabidopsis.h5.bz2


USAGE
Python PLEK2.py -i fasta_file -m model(ve: vertebrate , pl: plant)
#-fasta        The name of a fasta file, its sequences are to be predicted.
   
Examples:
$ python PLEK2.py -i PLEK2_test.fa -m ve

==============
Haotian Zhou, Master
School of Computer Science and Engineering,
Xi'an University of Technology,
5 South Jinhua Road,
Xi'an, Shaanxi 710048, P.R China

Aimin Li, PhD
School of Computer Science and Engineering,
Xi'an University of Technology,
5 South Jinhua Road,
Xi'an, Shaanxi 710048, P.R China

2350837044@qq.com 

liaiminmail AT gmail.com
emanlee815 AT 163.com
