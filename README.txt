The program recapSAM.py it’s a master degree project, with a free licence. This python program has to purpose to analyse SAM file (result of reads alignments on a reference sequence) and to offer a view of alignement quality to the user.

To use it you must execute recapSAM.py with on first argument the name of file to analyse (for exemple on a terminal : python3 recapSAM.py D:\System\Project\out.sam ).
The analyzed file can’t be a directory or a text file. That can’t be empty.
The rest of options (change the minimal value of read quality and export the result of analyze on a new text file) are proposed during the analyze.

RecapSAM.py give you som information on your SAM file :
-Number of reads 
-Rate of mapped reads
-Rate of mapped reads with a minimal quality (you can fix it)
-Number of pairs proper mapped reads
-Rate of totaly mapped reads with the minimal quality
