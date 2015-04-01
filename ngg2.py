#!/usr/bin/python

##########
# Import #
##########
import argparse
import re
import sys
import time

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split( '/' )[-1].split( '\\' )[-1]
VERSION = '1.0.0'

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
	
###########
# classes #
###########
class ContigIndex:
	def __init__( self, string ):
		vals = string.split( '\t' )
		
		self.name = vals[0]
		self.size = int( vals[1] )
		self.byteStart = int( vals[2] )
		self.basesPerLine = int( vals[3] )
		self.bytesPerLine = int( vals[4] )
		self.extraBytesPerLine = self.bytesPerLine - self.basesPerLine
	
	def getBytePosition( self, position ):
		if position < 1 or position > self.size:
			return None
		
		generalOffset = self.byteStart + position - 1
		sequenceLineWraps = int( ( position - 1 ) / self.basesPerLine )
		whiteSpaceChars = self.extraBytesPerLine * sequenceLineWraps
		
		return generalOffset + whiteSpaceChars

if __name__ == '__main__':
	############
	# argparse #
	############
	parser = argparse.ArgumentParser( prog=SCRIPT_NAME, epilog="%s v%s" % ( SCRIPT_NAME, VERSION ) )
	
	parser.add_argument( 'fastaFile', help="Path to indexed FASTA file to be scanned for gRNA sites" )
	parser.add_argument( '--outputFile', help="Defaults to gRNAs.csv", default="gRNAs.csv" )
	parser.add_argument( '--fastaIndex', help="Expects the index to be the full FASTA name plus '.fai'. A different name can be specified here.", default=None )
	parser.add_argument( '--region', help="Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig.", action='append' )
	parser.add_argument( '--allContigs', help="Ignores any defined regions and finds all gRNA sites in all contigs", default=False, action='store_true' )
	parser.add_argument( '--allowN', help="Allows N bases in gRNA site", action="store_true", default=False )
	parser.add_argument( '--onlyGstart', help="Only retain sites that start with G base", action="store_true", default=False )
	parser.add_argument( "--quiet", default=True, action='store_false', dest='verbose' )

	args = parser.parse_args()
	
	if args.fastaIndex == None:
		args.fastaIndex = "%s.fai" % ( args.fastaFile )
	
	if args.region is None and not args.allContigs:
		error_msg( "No regions were specified, and not set to check all contigs. Please specify a region or add '--allContigs' to your command" )
	elif args.allContigs:
		args.region = []

	options = "%s v%s\n\nOptions\n=======\n" % ( SCRIPT_NAME, VERSION )
	options += "FASTA: %s\n" % ( args.fastaFile )
	options += "FASTA Index: %s\n" % ( args.fastaIndex )
	options += "Output file: %s\n" % ( args.outputFile )
	options += "Region or all contigs: %s\n" % ( "Region" if not args.allContigs else "All contigs" )
	options += "Allow Ns in gRNA site: %s\n" % ( str( args.allowN ) )
	options += "Only G starts: %s\n" % ( str( args.onlyGstart ) )
	options += "Verbose: %s\n" % ( str( args.verbose ) )

	if args.verbose:
		log_msg( options )

	#################
	# Compile RegEx #
	#################
	whiteSpace = re.compile( '\n|\r' )
	
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
	
	########################
	# Read the FASTA index #
	########################
	contigIndexDict = {}
	
	try:
		FAIFH = open( args.fastaIndex, 'r' )
	except:
		error_msg( "Couldn't open FASTA index [%s]" % ( args.fastaIndex ) )
		
	for line in FAIFH:
		line = line.rstrip()
		
		if len( line ) == 0 or line[0] == '#':
			continue
		else:
			contigData = ContigIndex( line )
			
			if contigData.name in contigIndexDict:
				error_msg( "Contig [%s] listed twice in FASTA index!" % ( contigData.name ) )
			else:
				contigIndexDict[ contigData.name ] = contigData
	
	FAIFH.close()
	
	###################
	# Open FASTA file #
	###################
	try:
		FAFH = open( args.fastaFile, 'r' )
	except:
		error_msg( "Could not open FASTA file [%s]" % ( args.fastaFile ) )
	
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
	if args.allContigs:
		for name in contigIndexDict:
			args.region.append( '%s:1-%s' % ( name, contigIndexDict[name].size ) )
	
	for reg in args.region:
		vals = reg.split( ":" )
		
		if len( vals ) == 1:
			if vals[0] not in contigIndexDict:
				error_msg( "The contig [%s] was not indexed in the FASTA index file [%s]!" % ( vals[0], args.fastaIndex ) )
			
			contig = vals[0]
			start = 1
			end = contigIndexDict[contig].size
		else:
			contig = vals[0]
			pos = vals[1].split( '-' )
			start = int( pos[0] )
			end = int( pos[1] )
		
		if start > end:
			error_msg( "Start position [%s] is less than end position [%s]" % ( start, end ) )
		elif end - start + 1 < 23:
			error_msg( "Region too small to find gRNA (%s bp )" % ( end - start + 1 ) )
		elif contig not in contigIndexDict:
			error_msg( "Contig [%s] not found in FASTA index [%s]!" % ( contig, args.fastaIndex ) )
			
		#################
		# Pull sequence #
		#################
		startBytes = contigIndexDict[contig].getBytePosition( start )
		
		if startBytes == None:
			error_msg( "Start position was %s of contig [%s]" % ( 'before start' if start < 1 else 'after end', start ) )
		
		endBytes = contigIndexDict[contig].getBytePosition( end )
		
		if endBytes == None:
			error_msg( "End position was %s of contig [%s]" % ( 'before start' if end < 1 else 'after end', end ) )
		
		totalBytes = endBytes - startBytes + 1
		
		FAFH.seek( startBytes, 0 )
		
		sequence = FAFH.read( totalBytes )
		sequence = whiteSpace.sub( '', sequence ).upper()
		
		#################
		# Sense matches #
		#################
		for m in ngg2re.finditer( sequence ):
			grna = m.group()[:20]
			pam = m.group()[20:]
			startPos = start + m.start()
			endPos = startPos + len( grna ) - 1
			OUTFH.write( "%s,%s,%s,%s,%s,sense,%s\n" % ( contig, startPos, endPos, grna, pam, 'Yes' if grna[0] == 'G' else 'No' ) )
		
		#####################
		# Antisense matches #
		#####################
		for m in ccn2re.finditer( sequence ):
			grna = rev_comp( m.group()[3:] )
			pam = rev_comp( m.group()[:3] )
			startPos = start + m.start() + 3
			endPos = startPos + len( grna ) - 1
			OUTFH.write( "%s,%s,%s,%s,%s,antisense,%s\n" % ( contig, startPos, endPos, grna, pam, 'Yes' if grna[0] == 'G' else 'No' ) )
		
	OUTFH.close()
	FAFH.close()
