def main():
    fileIn1=open("14-Cyclonic-25m-B-R1-unmerged.fastq","r")
    fileOut1=open("14-Cyclonic-25m-B-R1-unmerged-removed.fastq","w")
    for i, line in enumerate(fileIn1):
        if i % 4 == 0:
            line = line.strip("\n")
            line0 = line
        elif i % 4 == 1:
            line = line.strip("\n")
            line1 = line
            length1 = len(line)
        elif i % 4 == 2:
            line = line.strip("\n")
            line2 = line
        elif i % 4 == 3:
            line = line.strip("\n")
            length2 = len(line)
            line3 = line
            if length1 == length2:
                print(line0, file=fileOut1)
                print(line1, file=fileOut1)
                print(line2, file=fileOut1)
                print(line3, file=fileOut1)
    fileIn1.close()
    fileOut1.close()


main()
