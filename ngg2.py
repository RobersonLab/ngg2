#!/usr/bin/python

##########
# import #
##########
import argparse
import multiprocessing as mp
import pyfaidx
import regex
import re
import sys
import time

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split( '/' )[-1].split( '\\' )[-1]
VERSION = 'development'

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

def compile_regex_patterns( type, allow_noncanonical ):
	if type == 'exhaustive':
		if allow_noncanonical == False:
			sense_regex = regex.compile( 'G[ACGT]{17}GG[ACGT]{1}GG' )
			antisense_regex = regex.compile( 'CC[ACGT]{1}CC[ACGT]{17}C' )
		elif allow_noncanonical == True:
			sense_regex = regex.compile( '[ACGT]{18}GG[ACGT]{1}GG' )
			antisense_regex = regex.compile( 'CC[ACGT]{1}CC[ACGT]{18}' )
		else:
			error_msg( "Only logical accepted in compile_regex_patterns when describing only G starts" )
	
	elif type == 'block':
		if allow_noncanonical == False:
			sense_regex = re.compile( 'G[ACGT]{17}GG[ACGT]{1}GG' )
			antisense_regex = re.compile( 'CC[ACGT]{1}CC[ACGT]{17}C' )
		elif allow_noncanonical == True:
			sense_regex = re.compile( '[ACGT]{18}GG[ACGT]{1}GG' )
			antisense_regex = re.compile( 'CC[ACGT]{1}CC[ACGT]{18}' )
		else:
			error_msg( "Only logical accepted in compile_regex_patterns when describing only G starts" )
	
	else:
		error_msg( "Type [%s] not recognized." % ( type ) )
	
	return ( sense_regex, antisense_regex )

def tuple_to_key( tuple ):
	return '_'.join( [ str(x) for x in tuple ] )

###########
# classes #
###########
class Grna:
	"""The positionStart is expected to be 0-base. Very important!"""
	def __init__( self, contig, positionStart, strand, regexMatch ):
		if strand == '+':
			self.seq = regexMatch.group()[:20]
			self.pam = regexMatch.group()[20:]
			self.start = positionStart + regexMatch.start() + 1 # switch to 1-base
			self.end = self.start + len( self.seq ) - 1
		elif strand == '-':
			self.seq = rev_comp( regexMatch.group()[3:] )
			self.pam = rev_comp( regexMatch.group()[:3] )
			self.start = positionStart + regexMatch.start() + 3 + 1 # switch to 1-base
			self.end = self.start + len( self.seq ) - 1
		else:
			error_msg( "Strand [%s] not recognized! Please specify '+' or '-'" )
		
		self.gCount = self.seq.count( 'G' )
		self.contig = contig
		self.strand = strand
		self.gStart = True if self.seq[0] == 'G' else False
	
	def __str__( self ):
		return "%s,%s,%s,%s,%s,%s,%s" % ( self.contig, self.start, self.end, self.seq, self.pam, self.strand, self.gStart )
		   
	def __repr__( self ):
		return str( self )

######################
# fxns using classes #
######################
def unbuffered_ngg_scan( fasta_fname, output_fname, region_list, regex_tuple, scan_type, max_g_content ):
	"""The locusTuple is expected in 0-start begin and 1-start end from pyfaidx. If you generate the tuple yourself, make sure it is in that format."""
	
	ngg2re, ccn2re = regex_tuple
	count = 0
	
	# the list of regions includes plus and neg strand for multiproc method
	# so we need to just keep half of them
	region_list = [region_list[x] for x in range(0, len(region_list),2)]
	
	with pyfaidx.Fasta( fasta_fname, as_raw=True ) as FAIDX, open( output_fname, 'w' ) as OUTFH:
		OUTFH.write( "Contig,Start,End,gRNA_Seq,PAM,Strand,G_start\n" )
		
		for currTuple in region_list:
			contig, start, end, strand = currTuple
		
			sequence = str( FAIDX[contig][start:end] ).upper()
			
			#################
			# Sense matches #
			#################
			if scan_type == 'block':
				for m in ngg2re.finditer( sequence ):
					grna = Grna( contig, start, '+', m )
					
					if grna.gCount <= max_g_content:
						count += 1
						OUTFH.write( "%s\n" % ( grna ) )
			else:
				for m in ngg2re.finditer( sequence, overlapped=True ):
					grna = Grna( contig, start, '+', m )
					
					if grna.gCount <= max_g_content:
						count += 1
						OUTFH.write( "%s\n" % ( grna ) )
					
			#####################
			# Antisense matches #
			#####################
			if scan_type == 'block':
				for m in ccn2re.finditer( sequence ):
					grna = Grna( contig, start, '-', m )
					
					if grna.gCount <= max_g_content:
						count += 1
						OUTFH.write( "%s\n" % ( grna ) )
			else:
				for m in ccn2re.finditer( sequence, overlapped=True ):
					grna = Grna( contig, start, '-', m )
					
					if grna.gCount <= max_g_content:
						count += 1
						OUTFH.write( "%s\n" % ( grna ) )

	return count

def multiproc_ngg_scan( fasta_fname, region_tuple, regex_tuple, scan_type, max_g_content ):
	"""Again the region_tuple start is expected to be 0-based, and end is expected to be 1-based"""
	
	ngg2re, ccn2re = regex_tuple
	
	contig, start, end, strand = region_tuple
	
	site_list = []
	
	with pyfaidx.Fasta( fasta_fname, as_raw=True ) as FAIDX:
		sequence = str( FAIDX[contig][start:end] ).upper()
	
	if strand == '+':
		#################
		# Sense matches #
		#################
		if scan_type == 'block':
			for m in ngg2re.finditer( sequence ):
				grna = Grna( contig, start, '+', m )
				
				if grna.gCount <= max_g_content:
					site_list.append( grna )
		else:
			for m in ngg2re.finditer( sequence, overlapped=True ):
				grna = Grna( contig, start, '+', m )
				
				if grna.gCount <= max_g_content:
					site_list.append( grna )
	
	elif strand == '-':
		#####################
		# Antisense matches #
		#####################
		if scan_type == 'block':
			for m in ccn2re.finditer( sequence ):
				grna = Grna( contig, start, '-', m )
				
				if grna.gCount <= max_g_content:
					site_list.append( grna )
		else:
			for m in ccn2re.finditer( sequence, overlapped=True ):
				grna = Grna( contig, start, '-', m )
				
				if grna.gCount <= max_g_content:
					site_list.append( grna )
	else:
		error_msg( "Strand [%s] not recognized" % ( strand ) )
	
	return ( tuple_to_key( region_tuple ), site_list )
	
########
# main #
########
if __name__ == '__main__':
	############
	# argparse #
	############
	parser = argparse.ArgumentParser( prog=SCRIPT_NAME, epilog="%s v%s" % ( SCRIPT_NAME, VERSION ) )

	# FASTA index will be created if it does not exist when pyfaidx Fasta is initialized
	parser.add_argument( 'fastaFile', help="Path to FASTA file to be scanned for gRNA sites" )
	parser.add_argument( '--outputFile', help="Defaults to gRNAs.csv", default="gRNAs.csv" )
	parser.add_argument( '--region', help="Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire genome will be searched!", action='append', required=False )
	parser.add_argument( '--allowNoncanonical', help="Allow sites that start with bases other than G", action="store_true", default=False )
	parser.add_argument( '--blockScan', help="Fast, but imperfect scanning that ignores gRNA sites that does not find overlapping gRNA sites", action='store_true', default=False )
	parser.add_argument( '--skipUniqueScan', help="Do not test for uniqueness if invoked (can still use multiple cores)", default=False, action="store_true" )
	parser.add_argument( '--unbuffered', help="Scan contigs in order with one processor and immediately write results. Block scan + unbuffered may provide fasta performance", default=False, action="store_true" )
	parser.add_argument( '--onlyUnique', help="Only report unique hits in output and counts", default=False, action="store_true" )
	parser.add_argument( '--maxSiteGs', help="ATTN", default=15, type=int )
	parser.add_argument( '--cores', help="Set to run contigs / strands on multiple cores simultaneously", type=int, default=1 )
	parser.add_argument( "--quiet", default=True, action='store_false', dest='verbose' )

	args = parser.parse_args()
	
	##############################
	# fix up conflicting options #
	##############################
	if args.unbuffered:
		args.cores = 1
		args.skipUniqueScan = True
		args.onlyUnique = False
	elif args.skipUniqueScan:
		args.onlyUnique = False
	
	#### ATTN SET REGION VS ALL
	
	######################
	# set options string #
	######################
	options = "%s v%s\n\nOptions\n=======\n" % ( SCRIPT_NAME, VERSION )
	options += "FASTA: %s\n" % ( args.fastaFile )
	options += "Output file: %s\n" % ( args.outputFile )
	options += "Target: %s\n" % ( "Region(s)" if not args.region else "All contigs" )
	options += "Allow non-canonical starts?: %s\n" % ( str( args.allowNoncanonical ) )
	options += "Max G-bases per site: %s\n" % ( args.maxSiteGs )# max G content
	options += "Scan type: %s\n" % ( "Block" if args.blockScan else "Exhaustive" ) # scan type
	options += "Buffered scan: %s\n" % ( "No" if args.unbuffered else "Yes" ) # unbuffered
	options += "Test site uniqueness: %s\n" % ( "No" if args.unbuffered or args.skipUniqueScan else "Yes" ) # skip unique scan
	options += "Only unique sites: %s\n" % ( "Yes" if not args.unbuffered and not args.skipUniqueScan and args.onlyUnique else "No" ) # only unique
	options += "Processes: %s\n" % ( args.cores ) # cores
	options += "Verbose: %s\n" % ( str( args.verbose ) )
	
	if args.verbose:
		log_msg( options )

	###################
	# Open FASTA file #
	###################
	region_list = []
	
	with pyfaidx.Fasta( args.fastaFile, as_raw=True ) as FAIDX:
		if args.region is None:
			args.region = FAIDX.keys()
		for reg in args.region:
			contig, start, end = pyfaidx.ucsc_split( reg )  # returns [0,1) coordinates with NoneType if not start or end

			# no need for checking for contig in FASTA file, as this is done by pyfaidx
			if start is None:
				start = 1
			if end is None:
				end = len( FAIDX[contig] )

			if end > len( FAIDX[contig] ):
				error_msg( "End position [%s] is past the end of [%s] (%s)" % ( end, contig, len( FAIDX[contig] ) ) )
			elif start > end:
				error_msg( "Start position [%s] is less than end position [%s]" % ( start, end ) )
			elif end - start < 23:  # [0, 1) coordinates make arithmetic easier
				error_msg( "Region too small to find gRNA (%s bp )" % ( end - start ) )
			if contig not in FAIDX:
				error_msg( "Contig [%s] not found in FASTA index [%s]!" % ( contig, FAIDX.filename ) )
				
			region_list.append( ( contig, start, end, '+' ) )
			region_list.append( ( contig, start, end, '-' ) )
	
	#################
	# Compile RegEx #
	#################
	regexTuple = compile_regex_patterns( 'block' if args.blockScan else 'exhaustive', args.allowNoncanonical )

	################
	# process data #
	################
	site_count = 0
	unique_count = 0
	
	if args.unbuffered:
		site_count = unbuffered_ngg_scan( args.fastaFile, args.outputFile, region_list, regexTuple, "block" if args.blockScan else 'exhaustive', args.maxSiteGs )
	else:
		#############################
		# async process the contigs #
		#############################
		work_pool = mp.Pool( processes=args.cores )
		results = [ work_pool.apply_async( multiproc_ngg_scan, args=(args.fastaFile, x, regexTuple, "block" if args.blockScan else 'exact', args.maxSiteGs ) ) for x in region_list ]
		grna_output = [ x.get() for x in results ]
		
		####################
		# figure order out #
		####################
		result_order = {}
		for index in range( len( grna_output ) ):
			result_order[grna_output[index][0]] = index
				
		##########################
		# write outputs in order #
		##########################
		with open( args.outputFile, 'w' ) as OUTFH:
			if args.skipUniqueScan:
				OUTFH.write( "Contig,Start,End,gRNA_Seq,PAM,Strand,G_start\n" )
				
				for region in region_list:
					region_key = tuple_to_key( region )
					region_index = result_order[region_key]
					
					site_count += len( grna_output[region_index][1] )
				
					for grna in grna_output[region_index][1]:
						OUTFH.write( "%s\n" % ( grna ) )
				
			else:
				################################
				# dictionary of sites to count #
				################################
				grna_dictionary = {}
				get = grna_dictionary.get
				
				for curr_grna_tuple in grna_output:
					for curr_grna in curr_grna_tuple[1]:
						grna_dictionary[curr_grna.seq] = get( curr_grna.seq, 0 ) + 1
				
				################
				# Write header #
				################
				OUTFH.write( "Contig,Start,End,gRNA_Seq,PAM,Strand,G_start,Unique\n" )
				
				if args.onlyUnique:
					###########################
					# Only print unique gRNAs #
					###########################
					for region in region_list:
						region_key = tuple_to_key( region )
						region_index = result_order[region_key]
						
						site_count += len( grna_output[region_index][1] )
					
						for grna in grna_output[region_index][1]:
							if grna_dictionary[grna.seq] == 1:
								unique_count += 1
								OUTFH.write( "%s,Yes\n" % ( grna ) )
				else:
					############################################
					# or print everything, noting unique gRNAs #
					############################################
					for region in region_list:
						region_key = tuple_to_key( region )
						region_index = result_order[region_key]
						
						site_count += len( grna_output[region_index][1] )
					
						for grna in grna_output[region_index][1]:
							if grna_dictionary[grna.seq] == 1:
								unique_count += 1
								OUTFH.write( "%s,Yes\n" % ( grna ) )
							else:
								OUTFH.write( "%s,No\n" % ( grna ) )
	
	log_msg( "%s total sites found by %s scanning" % ( site_count, "block" if args.blockScan else "exhaustive" ), printTime = False )
	if not args.skipUniqueScan and not args.unbuffered:
		log_msg( "%s unique sites (%.1f%%)" % ( unique_count, float( unique_count ) / float( site_count ) * 100.0 ), printTime = False )
