import pysam
import os
import sys


def filter_bam(in_bam,out_bam1,out_bam2):
    """
    """

    IN = pysam.AlignmentFile(in_bam,"rb")
    OUT1 = pysam.AlignmentFile(out_bam1,"wb",template=IN)
    OUT2 = pysam.AlignmentFile(out_bam2,"wb",template=IN)
    count1=0
    count2=0
    
    for read in IN:
        tags = read.tags
        for tag in tags:
            if tag[0] == "cD": ## Consensus depth
                if tag[1] > 1:
                    OUT1.write(read)
                    count1+=1
                if tag[1] > 2:
                    OUT2.write(read)
                    count2+=1
    IN.close()
    OUT1.close()
    OUT2.close()
    print "Filtered : \n min_reads=2 : %i\n min_reads=3 : %i\n"%(count1,count2)
    
if __name__ == "__main__":
    filter_bam(sys.argv[1],sys.argv[2],sys.argv[3])
