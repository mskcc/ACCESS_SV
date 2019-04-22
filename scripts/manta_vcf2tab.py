"""
dellyVcf2Tab
~~~~~~~~~~~~

:Description: This module converts the Delly Vcf file having tumor normal, to tab-delimited format for input to iAnnotateSV

"""
'''
Created on Mar 18, 2015
Description: This module converts the Delly Vcf file having tumor normal, to tab-delimited format for input to iAnnotateSV
@author: Ronak H Shah
::Input::
vcfFile: Input vcf file to convert
outputFileName: Name of the output file
OutputDir: Directory for output file
::Output::
outputFile: Tab-delimited file containing:
chr1: Its the chromosome name for first break point [1,2,3,4,5,6,7 etc..],
pos1: Its the chromosome loaction for first break point [1-based],
str1: Its the read direction for the first break point [0=top/plus/reference, 1=bottom/minus/complement],
chr2: Its the chromosome name for second break point [1,2,3,4,5,6,7 etc..],
pos2: Its the chromosome loaction for second break point [1-based],
str2: Its the read direction for the second break point [0=top/plus/reference, 1=bottom/minus/complement],
'''

import os
import vcf
# import checkparameters as cp
import logging
import coloredlogs
import argparse



def vcf2tab(vcfFile, outputDir, verbose):
    logger = logging.getLogger('iCallSV.dellyVcf2Tab')
    coloredlogs.install(level='DEBUG')
    """This ``converts`` the Delly Vcf file having tumor normal, to tab-delimited format for input to iAnnotateSV


    :param str vcfFile: str of vcf file to be converted
    :param str outputDir: str for the output directory
    :param bool verbose: a boolean
    :return: A str name of tab-delimited file
    :rtype: str

    """
    # cp.checkFile(vcfFile)
    # cp.checkDir(outputDir)
    if(verbose):
        logger.info("dellyVcf2Tab: All Input Parameters look good. Lets convert to tab-delimited file")
    vcf_reader = vcf.Reader(open(vcfFile, 'r'))
    outputFileName = os.path.splitext((os.path.basename(vcfFile)))[0] + ".tab"
    outputFile = outputDir + "/" + outputFileName
    outputHandle = open(outputFile, "w")
    outputHandle.write("chr1\tpos1\tstr1\tchr2\tpos2\tstr2\n")
    for record in vcf_reader:
        (chrom1,
         start1,
         start2,
         chrom2,
         contype,
         str1,
         str2) = (None for i in range(7))
        chrom1 = record.CHROM
        start1 = record.POS
        if("END" in record.INFO):
            start2 = record.INFO['END']
        if("CHR2" in record.INFO):
            chrom2 = str(record.INFO['CHR2'][0])
        if("CT" in record.INFO):
            contype = str(record.INFO['CT'][0])
        # print contype
        (startCT, endCT) = contype.split("to")
        if((int(startCT) == 3) and (int(endCT) == 3)):
            str1 = 0
            str2 = 0
        elif((int(startCT) == 3) and (int(endCT) == 5)):
            str1 = 0
            str2 = 1
        elif((int(startCT) == 5) and (int(endCT) == 3)):
            str1 = 1
            str2 = 0
        elif((int(startCT) == 5) and (int(endCT) == 5)):
            str1 = 1
            str2 = 1
        else:
            if(verbose):
                logger.info(
                    "dellyVcf2Tab: The connection type (CT) given in the vcf file is incorrect.CT: %s",
                    contype)
        outputHandle.write(
            str(chrom1) +
            "\t" +
            str(start1) +
            "\t" +
            str(str1) +
            "\t" +
            str(chrom2) +
            "\t" +
            str(start2) +
            "\t" +
            str(str2) +
            "\n")
    outputHandle.close()
    if(verbose):
        logger.info("dellyVcf2Tab: Finished conversion of Vcf file to tab-delimited file")
        logger.info("dellyVcf2Tab: Output can be found: %s", outputFile)
    return(outputFile)

if __name__ == '__main__':
    # Setup argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="set verbosity level [default: %(default)s]")
    parser.add_argument(
        "-i",
        "--inputvcf",
        action="store",
        dest="input_vcf",
        required=True,
        metavar='/some/path/to/vcf.vcf',
        help="Full path to the vcf file.")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outdir",
        required=True,
        metavar='/somepath/output',
        help="Full Path to the output dir.")

    args = parser.parse_args()
    vcf2tab(args.input_vcf, args.outdir, args.verbose)
