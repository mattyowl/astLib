#!/usr/bin/python
#File: coordConvert/coordConvert.py
#Created: Sat Dec 15 17:08:26 2012
#Last Change: Sat Dec 15 17:15:07 2012
# -*- coding: utf-8 -*-
#
# Converts file fed in into hms, dms format
#

from astLib import astCoords
import sys

# Main
if len(sys.argv) < 8:
    print ("Run: % coordConvert.py <in file> <out file> <RA column> <dec column> <way [either \"to_decimal\" or \"to_hmsdms\"]> <coord delimiter> <output column delimiter [e.g. \"tab\" \",\"\"&\" etc.]>")
else:

    inFile = sys.argv[1]
    outFile = sys.argv[2]
    RACol = int(sys.argv[3])-1
    decCol = int(sys.argv[4])-1
    way = sys.argv[5]
    delimiter = sys.argv[6]
    colDelim = str(sys.argv[7]).rstrip("\n")

    if colDelim == "tab":
        colDelim = "\t"

    wayOk = False
    if way == "to_decimal":
        wayOk = True
    if way == "to_hmsdms":
        wayOk = True
    if wayOk == False:
        print("<way>: either \"to_decimal\" or \"to_hmsdms\"")
        sys.exit()

    reader = open(inFile, "r")
    lines = reader.readlines()
    reader.close()

    writer = open(outFile, "w")
    for row in lines:
        if len(row)>1:
            if "#" not in row[:1]:
                rowBits = row.split("\t")
                if way == "to_decimal":
                    RADeg = astCoords.hms2decimal(rowBits[RACol], delimiter)
                    decDeg = astCoords.dms2decimal(rowBits[decCol], delimiter)
                if way == "to_hmsdms":
                    RADeg = astCoords.decimal2hms(float(rowBits[RACol]),
                            delimiter)
                    decDeg = astCoords.decimal2dms(float(rowBits[decCol]),
                            delimiter)
                writeString = ""
                for i in range(len(rowBits)):
                    if i == RACol:
                        writeString = writeString+str(RADeg)+"\t"
                    elif i == decCol:
                        writeString = writeString+str(decDeg)+"\t"
                    elif rowBits[i].find("\n") != -1:
                        writeString = writeString+str(rowBits[i])
                    else:
                        writeString = writeString+str(rowBits[i])+"\t"
                # new line character already included
                writer.write(writeString.replace("\n", "").replace("\t",
                    colDelim)+"\n")
            else:
                writer.write(row.replace("\t", colDelim))
    writer.close()

#-----------------------------------------------------------------------------
