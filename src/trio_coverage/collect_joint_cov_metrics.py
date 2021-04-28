'''
Name: collect_joint_cov_metrics.py
Authors: Matt Halvorsen, adapted from code written by Josh Bridgers
Desc: Create a summary matrix of called real estate across a population
      of samples. Generate summary statistics across samples and across
      CCDS genes.
'''

import os
import sys
import argparse
import datetime
import glob
import math
import gzip

from scipy import stats

def main():
    ''' main function for analysis of coverage in trio and output of cov stat matrix '''
    desc = "Create a summary matrix of called real estate across a population of\
    samples.  Generate summary statistics and inherent QC steps."

    usage = "python2.7 collect_joint_cov_metrics.py " + \
            "-r [revision number] <trios.covpath> <mat_out_filename>"

    parser = argparse.ArgumentParser(prog="collect_joint_cov_metrics", 
                                     usage=usage, description=desc)

    parser.add_argument('-r', action="store", dest="revVer", type=int, 
                        default='14', choices=[9,14,20], 
                        help='CCDS Revision Version [9,14,20]')
    parser.add_argument('covpath_filename', action="store", 
                        help='file with list of IGM joint cov files')
    parser.add_argument('mat_filename', action="store", 
                        help='name of trio cov matrix file')
    args = parser.parse_args()

    # checkFnames(args.covpath_filename, args.revVer)
    createMatrix(args.covpath_filename, args.mat_filename)
    return

def openFile(filename):
    '''Open filename in a way that minds filetype'''
    if filename.find(".gz") != -1:
        fh = gzip.open(filename, "rb")
    else:
        fh = open(filename, "r")
    return fh

def createMatrix(filenames,mat_filename):
    '''Creates a gene matrix by parsing information from the list of files'''
    fileCount=0
    geneDict = {}
    sampleDict = {}
    allFileNames = []

    for file in open(filenames).readlines():
        sampleName =file.strip().split('/')[-1].split('.')[0]
        allFileNames.append(sampleName)
        fileCount+=1

        totalGeneSize = 0
        coveredGeneSize = 0
        excYTotalGeneSize = 0
        excYCoveredGeneSize = 0

        if fileCount == 1:
            fh = openFile(file.strip())
            for genes in fh.readlines():
                geneInfo = genes.strip().split()

                ## NEW 
                if int(geneInfo[3]) > int(geneInfo[4]):
                    geneInfo[3] = geneInfo[4]
                
                if geneInfo[0] != 'X' and geneInfo[0] != 'Y':
                    geneDict[geneInfo[0].zfill(2)+'_'+geneInfo[1]] = [geneInfo[0],geneInfo[1],geneInfo[4],geneInfo[3]]
                else:
                    geneDict[geneInfo[0]+'_'+geneInfo[1]] = [geneInfo[0],geneInfo[1],geneInfo[4],geneInfo[3]]

                totalGeneSize+=int(geneInfo[4])
                coveredGeneSize+=int(geneInfo[3])

                if geneInfo[0] != 'Y':
                    excYTotalGeneSize+=int(geneInfo[4])
                    excYCoveredGeneSize+=int(geneInfo[3])

            geneCount = len(geneDict.keys())
            #print sampleName,coveredGeneSize,totalGeneSize,excYTotalGeneSize,excYCoveredGeneSize

            sampleDict[sampleName] = [coveredGeneSize,totalGeneSize,excYCoveredGeneSize,excYTotalGeneSize]

        else:

            otherGeneCount = 0
            fh = openFile(file.rstrip())
            for genes in fh.readlines():
                otherGeneCount += 1
                geneInfo = genes.strip().split()

                ## NEW 
                if int(geneInfo[3]) > int(geneInfo[4]):
                    geneInfo[3] = geneInfo[4]
                
                if geneInfo[0] != 'X' and geneInfo[0] != 'Y':
                    geneDict[geneInfo[0].zfill(2)+'_'+geneInfo[1]].append(geneInfo[3])
                else:
                    geneDict[geneInfo[0]+'_'+geneInfo[1]].append(geneInfo[3])

                totalGeneSize+=int(geneInfo[4])
                coveredGeneSize+=int(geneInfo[3])

                if geneInfo[0] != 'Y':
                    excYTotalGeneSize+=int(geneInfo[4])
                    excYCoveredGeneSize+=int(geneInfo[3])

            #print sampleName,coveredGeneSize,totalGeneSize,excYTotalGeneSize,excYCoveredGeneSize
            sampleDict[sampleName] = [coveredGeneSize,totalGeneSize,excYCoveredGeneSize,excYTotalGeneSize]

            if otherGeneCount != geneCount:
                print geneCount,otherGeneCount
                raise Exception, "%s does not have the same number of genes.  File likely truncated" % file.strip()


    #print geneDict['11_SERPINH1']
    outFile = open(mat_filename+'.tsv','w')
    outFile.write('Chr\tHGNC/CCDS\tOriginalSize\t%s\n' % "\t".join(allFileNames))
    outFile2 = open(mat_filename+'.genicProfileSummary.tsv','w')
    outFile2.write('Chr\tHGNC/CCDS\tOriginalSize\tMean\tMedian\tMin\tMax\tStd Dev\t%Mean\t%Median\t%Min\t%Max\n')
    for key in sorted(geneDict.keys()):
        if len(geneDict[key]) != int(fileCount)+3:
            raise Exception, "Gene %s does not have the correct number of entries" % key
        else:
            outFile.write('\t'.join(geneDict[key])+'\n')
            results = sorted(map(int,geneDict[key][3:]))
            origSize = int(geneDict[key][2])
            mean = round(float(sum(results))/len(results),7)
            geneMin = min(results)
            geneMax = max(results)
            stdev = round(math.sqrt(sum(map(lambda x: (x-mean)**2,results))/float(len(results))),7)
            median = getMedian(results)
            perMean = round(mean/origSize * 100,8)
            perMedian = round(float(median)/origSize * 100,8)
            perMin = round(float(geneMin)/origSize * 100,8)
            perMax = round(float(geneMax)/origSize * 100,8)

            if mean == 0.0:
                mean = int(mean)
            if stdev == 0.0:
                stdev = int(stdev)
            if perMean == 0.0:
                perMean = int(perMean)
            if perMedian == 0.0:
                perMedian = int(perMedian)
            if perMin == 0.0:
                perMin = int(perMin)
            if perMax == 0.0:
                perMax = int(perMax)

            #q1 = results[int(math.ceil((len(results)*25)/100)-1)]
            #q2 = results[int(math.ceil((len(results)*50)/100)-1)]
            outFile2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ("\t".join(geneDict[key][0:3]),mean,median,geneMin,geneMax,stdev,perMean,perMedian,perMin,perMax))

    outFile.close()
    outFile2.close()

    outFile3 = open(mat_filename+'.sampleProfileSummary.tsv','w')
    outFile3.write('Trio SampleID\tCovered Bases (Y exc)\tPercent Covered (Y exc)\n')
    for file in open(filenames).readlines():
        sampleName =file.strip().split('/')[-1].split('.')[0]
        sampleCoveredBases = sampleDict[sampleName][2]
        sampleTotalBases = float(sampleDict[sampleName][3])
        outFile3.write('%s\t%s\t%s\n' % (sampleName,sampleCoveredBases,round(sampleCoveredBases/sampleTotalBases,9)))
    outFile3.close()
    return

def getMedian(results):
    ''' get median value from a user-provided list of numeric values '''
    sortedLst = sorted(results)
    lstLen = len(results)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def checkFnames(filenames,revVer):
    '''Sanity checks on the list of file names'''
    if len(glob.glob(filenames)) == 0:
        raise Exception, "List of files names not found!"

    for file in open(filenames).readlines():
        file2 = file.strip().split('.')
        if file2[-2] != "r" + str(revVer) and file2[-3] != "r" + str(revVer):
            raise Exception, "File %s is not the specified revision number!" % file.strip()
        if len(glob.glob(file.strip())) == 0:
            raise Exception, "File %s is not found!" % file.strip()
    return

if __name__ == '__main__':
    main()
