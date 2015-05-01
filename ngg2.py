#!/usr/bin/python

##########
# Import #
##########
import argparse
import re
import sys
import time
import pyfaidx

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split( '/' )[-1].split( '\\' )[-1]
VERSION = '1.1.0'

########
# fxns #
########
def log_msg( msg, printTime=True ):
    if printTime:
        print time.strftime( '%H:%M %z on %b %d, %Y' )
    print msg
    sys.stdout.flush()

def error_msg( msg ):
    log_msg( "Error: %s\n" % ( msg.rstrip() ) )
    sys.exit( 1 )

def rev_comp( seq ):
    nuc_dict = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N' }

    return ''.join( [ nuc_dict[c] for c in seq.upper()[::-1] ] )

if __name__ == '__main__':
    ############
    # argparse #
    ############
    parser = argparse.ArgumentParser( prog=SCRIPT_NAME, epilog="%s v%s" % ( SCRIPT_NAME, VERSION ) )

    # FASTA index will be created if it does not exist when pyfaidx.Fasta is initialized
    parser.add_argument( 'fastaFile', help="Path to FASTA file to be scanned for gRNA sites" )
    parser.add_argument( '--outputFile', help="Defaults to gRNAs.csv", default="gRNAs.csv" )
    parser.add_argument( '--region', help="Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire genome will be searched!", action='append' )
    parser.add_argument( '--allowN', help="Allows N bases in gRNA site", action="store_true", default=False )
    parser.add_argument( '--onlyGstart', help="Only retain sites that start with G base", action="store_true", default=False )
    parser.add_argument( "--quiet", default=True, action='store_false', dest='verbose' )

    args = parser.parse_args()

    ###################
    # Open FASTA file #
    ###################
    try:
        FAIDX = pyfaidx.Fasta( args.fastaFile, as_raw=True )
    except:
        error_msg( "Could not open FASTA file [%s]" % ( args.fastaFile ) )

    if args.region is None:
        args.region = FAIDX.keys()

    options = "%s v%s\n\nOptions\n=======\n" % ( SCRIPT_NAME, VERSION )
    options += "FASTA: %s\n" % ( args.fastaFile )
    options += "Output file: %s\n" % ( args.outputFile )
    options += "Region or all contigs: %s\n" % ( "Region" if not args.region else "All contigs" )
    options += "Allow Ns in gRNA site: %s\n" % ( str( args.allowN ) )
    options += "Only G starts: %s\n" % ( str( args.onlyGstart ) )
    options += "Verbose: %s\n" % ( str( args.verbose ) )

    if args.verbose:
        log_msg( options )

    #################
    # Compile RegEx #
    #################
    if args.allowN:
        if args.onlyGstart:
            ngg2re = re.compile( 'G[ACGTN]{17}GG[NACGT]{1}GG' )
            ccn2re = re.compile( 'CC[ACGTN]{1}CC[NACGT]{17}C' )
        else:
            ngg2re = re.compile( '[ACGTN]{18}GG[NACGT]{1}GG' )
            ccn2re = re.compile( 'CC[ACGTN]{1}CC[NACGT]{18}' )
    else:
        if args.onlyGstart:
            ngg2re = re.compile( 'G[ACGT]{17}GG[ACGT]{1}GG' )
            ccn2re = re.compile( 'CC[ACGT]{1}CC[ACGT]{17}C' )
        else:
            ngg2re = re.compile( '[ACGT]{18}GG[ACGT]{1}GG' )
            ccn2re = re.compile( 'CC[ACGT]{1}CC[ACGT]{18}' )

    ###############
    # Open output #
    ###############
    try:
        OUTFH = open( args.outputFile, 'w' )
    except:
        error_msg( "Could not open output file [%s] for writing" % ( args.outputFile ) )

    OUTFH.write( "Contig,Start,End,gRNA_Seq,PAM,Strand,G_start\n" )

    ##################################
    # Operate over specified regions #
    ##################################

    for reg in args.region:
        contig, start, end = pyfaidx.ucsc_split(reg)  # returns [0,1) coordinates with NoneType if not start or end

        # no need for checking for contig in FASTA file, as this is done by pyfaidx

        if start is None:
            start = 1
        if end is None:
            end = len(FAIDX[contig])

        if end > len(FAIDX[contig]):
            error_msg( "End position [%s] is past the end of [%s] (%s)" % ( end, contig, len(FAIDX[contig]) ) )
        elif start > end:
            error_msg( "Start position [%s] is less than end position [%s]" % ( start, end ) )
        elif end - start < 23:  # [0, 1) coordinates make arithmetic easier
            error_msg( "Region too small to find gRNA (%s bp )" % ( end - start ) )
        if contig not in FAIDX:
            error_msg( "Contig [%s] not found in FASTA index [%s]!" % ( contig, FAIDX.filename ) )

        #################
        # Pull sequence #
        #################

        sequence = str( FAIDX[contig][start:end] ).upper()

        #################
        # Sense matches #
        #################
        for m in ngg2re.finditer( sequence ):
            grna = m.group()[:20]
            pam = m.group()[20:]
            startPos = start + m.start() + 1 # output is 1-base, not 0-base bed
            endPos = startPos + len( grna ) - 1
            OUTFH.write( "%s,%s,%s,%s,%s,sense,%s\n" % ( contig, startPos, endPos, grna, pam, 'Yes' if grna[0] == 'G' else 'No' ) )

        #####################
        # Antisense matches #
        #####################
        for m in ccn2re.finditer( sequence ):
            grna = rev_comp( m.group()[3:] )
            pam = rev_comp( m.group()[:3] )
            startPos = start + m.start() + 3 + 1 # output is 1-base, not 0-base bed
            endPos = startPos + len( grna ) - 1
            OUTFH.write( "%s,%s,%s,%s,%s,antisense,%s\n" % ( contig, startPos, endPos, grna, pam, 'Yes' if grna[0] == 'G' else 'No' ) )

    OUTFH.close()
    FAIDX.close()
