#!/opt/conda/bin/python

# import libs
import csv
import sys, getopt
import os

def read_arguments(argv):

    insheet = ''
    outsheet = ''
    indextype = ''

    usage="> Usage: ctg-parse-samplesheet.10x.py -s INPUT-SHEET -o OUTPUT-SHEET -i INDEX-TYPE [ -h HELP ] \n\n> Required columns (with the following header names): \n [Lane,Sample_ID,index,Sample_Project]. \n - 'Lane' entries can be left blank if all lanes included. \n - 'index' entries must be ID (e.g. SI-TT-E4) if dual index."

    try:
        opts, args = getopt.getopt(argv,"hs:o:i:",["insheet=", "outsheet=", "indextype="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if len(sys.argv) <= 3:
        print("> Error: please specify all arguments:")
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-s", "--insheet"):
            insheet = arg
        elif opt in ("-o", "--outsheet"):
            outsheet = arg
        elif opt in ("-i", "--indextype"):
            indextype = arg   
            
    return insheet, outsheet, indextype

def main(argv):

    insheet, outsheet, index = read_arguments(argv)

    with open(outsheet, 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['[Data]'])
        
        if index == 'dual':
            writer.writerow(['Lane','Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project'])
        else:
            writer.writerow(['Lane','Sample_ID','index','Sample_Project'])
            
        with open(insheet, 'r') as infile:
            my_reader = csv.reader(infile, delimiter=',')
            # row counter to define first line
            row_idx=0
            for row in my_reader:
                # if first line - get index of the 3 columns needed
                if row_idx == 0:
                    laneidx = row.index('Lane')
                    sididx  = row.index('Sample_ID')
                    idxidx  = row.index('index')
                    projidx = row.index('Sample_Project')
                else:
                    currlane = row[laneidx]
                    currsid = row[sididx]
                    curridx = row[idxidx]
                    currproj = row[projidx]
                    
                        
                    if index == 'dual':
                        writer.writerow([currlane,currsid,currsid,'','',curridx,curridx,curridx,curridx,currproj])
                    else:
                        writer.writerow([currlane,currsid,curridx,currproj])

		   
                row_idx += 1


if __name__ == "__main__":
    main(sys.argv[1:])
