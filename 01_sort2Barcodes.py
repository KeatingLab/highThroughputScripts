#!~/anaconda3/envs/py2/bin/python

from Bio import SeqIO
import sys
import itertools 

def main():
    #Updated for corrected file paths by dw
    #Setup variables (could parse command line args instead)
    file_f = "/scratch/users/dswhit/DATA/titeseq_illumina/171220Kea_D17-10768_1_sequence.fastq"
    file_r = "/scratch/users/dswhit/DATA/titeseq_illumina/171220Kea_D17-10768_2_sequence.fastq"
#    file_f = "/scratch/users/dswhit/DATA/titeseq_illumina/test_1.fastq"
#    file_r = "/scratch/users/dswhit/DATA/titeseq_illumina/test_2.fastq"


    dir_out = "/scratch/users/dswhit/DATA/titeseq_illumina/workspace/"
    myFormat = "fastq-sanger"
    numBarcodes = (24*9) #1 is enrich sorts and controls, 2-9 is four TiteSeq runs
    
    records_f = SeqIO.parse(open(file_f,"rU"), myFormat)
    records_r = SeqIO.parse(open(file_r,"rU"), myFormat)
    myOpenOutStreams = openOutStreams(dir_out,numBarcodes)
    for (forward, reverse) in itertools.izip(records_f,records_r):   
#    for (forward, reverse) in zip(records_f,records_r):
        barcodeIndex = getBarcode(forward, reverse)
        if(barcodeIndex>-1):
            #print(barcodeIndex, forward.seq ,reverse.seq)
            SeqIO.write(forward, myOpenOutStreams[barcodeIndex], myFormat)
            SeqIO.write(reverse, myOpenOutStreams[barcodeIndex], myFormat)

    closeOutStreams(myOpenOutStreams)
    
    return 1
    
def openOutStreams(dir_out,num):
    myArray = []
    for i in range(num):
        print str(i)
        myArray.append(open(dir_out+"barcode_"+str(i),'w'))
    return myArray

def closeOutStreams(array):
    for each in array:
        each.close()
            
def getBarcode(forward, reverse):
    sequence = str(forward.seq)
    forward_Quality = forward.letter_annotations['phred_quality']
   
    #First check if barcode2 is available
    #Barcode2 = index (identifies which experiment)
    # Updated for 9 index sequences by dw
    myBarcodes2 = ["ATCACG","CGATGT","TTAGGC","TGACCA","ACAGTG","GCCAAT",
                    "CAGATC","ACTTGA","GATCAG"]
    myBarcodes2Index=-1
    
    #Trying to get the index from the rev complement of the start of the reverse read
    #Casting as a string removes the ability to use further features, obviously, try 
    #to do in reverse order
    # This worked, far as I can tell.
    rev_sequence = (reverse.seq)
    lastSix = rev_sequence[0:6]
    rv_lastSix = lastSix.reverse_complement()
    myBarcodes2Index = getClosestBarcode(rv_lastSix,myBarcodes2,5)  # 5 indicates score must be better than 5
    if(myBarcodes2Index==-1):
        return -1
    #Then check barcode1 quality
    for i in range(5):
        if forward_Quality[i]<20:
            return -1
    
    myBarcodes = ["ACTCG","ACTGT", "AATGC", "AGTCA", "ATACG", "ATAGC",
                    "CGATC", "CTAAG", "CTCGA", "CGAAT", "CTGGT", "CGGTT",
                    "GACTT", "GTTCA", "GATAC", "GAGCA", "GATGA", "GTCTG",
                    "TCGGA", "TGACC", "TACTG", "TCCAG", "TCGAC", "TAGCT"]
    firstFive = sequence[0:5]
    if(firstFive in myBarcodes):
        return (myBarcodes2Index*24)+myBarcodes.index(firstFive)
    else:
        return -1

def getClosestBarcode(seq,myBarcodes2,mustMatch):
    bestIndex = -1
    bestScore = -1
    mySeq = list(seq)
    
    for index in range(len(myBarcodes2)):
        mySum = 0
        for nt in range(6):
            if(mySeq[nt]==myBarcodes2[index][nt]):
                mySum+=1
        if(mySum>=bestScore):
            bestIndex=index
            bestScore=mySum
            
    if(bestScore>mustMatch):
        return bestIndex
    else:
        return -1

    
if __name__ == "__main__":
    sys.exit(main())
