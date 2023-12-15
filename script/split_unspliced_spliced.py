#!/maps/projects/dan1/apps/python/3.7.13/bin/python 
import sys 
import pysam

in_bam = sys.argv[1]
unspliced_bam = sys.argv[2]
spliced_bam = sys.argv[3] 

bam = pysam.AlignmentFile(in_bam, "rb")
unspliced= pysam.AlignmentFile(unspliced_bam, "wb", template=bam)
spliced = pysam.AlignmentFile(spliced_bam, "wb", template=bam)

for read in bam:
    if read.has_tag('XS'):
        spliced.write(read)
    else:
        unspliced.write(read)
    
unspliced.close()
spliced.close()
