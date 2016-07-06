#!/usr/bin/env python

##########
# import #
##########
import argparse
import logging
import multiprocessing as mp
import pyfaidx
import regex
import re
import sys
import time

from ngg2 import __script_path__, __script_name__, __version__, rev_comp, compile_regex_patterns, tuple_to_key, Grna, unbuffered_ngg_scan, multiproc_ngg_scan

def run():
	############
	# argparse #
	############
	parser = argparse.ArgumentParser( prog=__script_name__, epilog="%s v%s" % ( __script_name__, __version__ ) )

	# FASTA index will be created if it does not exist when pyfaidx Fasta is initialized
	parser.add_argument( 'fastaFile', help="Path to FASTA file to be scanned for gRNA sites" )
	parser.add_argument( '--outputFile', help="Defaults to gRNAs.csv", default="gRNAs.csv" )
	parser.add_argument( '--region', help="Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire genome will be searched!", action='append', required=False )
	parser.add_argument( '--allowNoncanonical', help="Allow sites that start with bases other than G", action="store_true", default=False )
	parser.add_argument( '--blockScan', help="Fast, but imperfect scanning that ignores gRNA sites that does not find overlapping gRNA sites", action='store_true', default=False )
	parser.add_argument( '--skipUniqueScan', help="Do not test for uniqueness if invoked (can still use multiple cores)", default=False, action="store_true" )
	parser.add_argument( '--unbuffered', help="Scan contigs in order with one processor and immediately write results. Block scan + unbuffered may provide fasta performance", default=False, action="store_true" )
	parser.add_argument( '--onlyUnique', help="Only report unique hits in output and counts", default=False, action="store_true" )
	parser.add_argument( '--maxSiteGs', help="Maximum # of G bases in gRNA site", default=15, type=int )
	parser.add_argument( '--cores', help="Set to run contigs / strands on multiple cores simultaneously", type=int, default=1 )
	parser.add_argument( "--loglevel", choices=[ 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL' ], default='INFO' )

	args = parser.parse_args()
	
	#################
	# setup logging #
	#################
	logging.basicConfig( format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
	logger = logging.getLogger( __script_name__ )
	logger.setLevel( args.loglevel )
	
	##############################
	# fix up conflicting options #
	##############################
	if args.unbuffered:
		args.cores = 1
		args.skipUniqueScan = True
		args.onlyUnique = False
	elif args.skipUniqueScan:
		args.onlyUnique = False
	
	######################
	# set options string #
	######################
	options = "%s v%s\n\nOptions\n=======\n" % ( __script_name__, __version__ )
	options += "FASTA: %s\n" % ( args.fastaFile )
	options += "Output file: %s\n" % ( args.outputFile )
	options += "Target: %s\n" % ( "All contigs" if args.region is None else str( args.region ) )
	options += "Allow non-canonical starts?: %s\n" % ( str( args.allowNoncanonical ) )
	options += "Max G-bases per site: %s\n" % ( args.maxSiteGs ) # max G content
	options += "Scan type: %s\n" % ( "Block" if args.blockScan else "Exhaustive" ) # scan type
	options += "Buffered scan: %s\n" % ( "No" if args.unbuffered else "Yes" ) # unbuffered
	options += "Test site uniqueness: %s\n" % ( "No" if args.unbuffered or args.skipUniqueScan else "Yes" ) # skip unique scan
	options += "Only unique sites: %s\n" % ( "Yes" if not args.unbuffered and not args.skipUniqueScan and args.onlyUnique else "No" ) # only unique
	options += "Processes: %s\n" % ( args.cores ) # cores
	
	logger.info( options )

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
				start = 0 # remember, pyfaidx expects 0 start
			if end is None:
				end = len( FAIDX[contig] )
				
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
	
	logger.info( "%s total sites found by %s scanning" % ( site_count, "block" if args.blockScan else "exhaustive" ) )
	if not args.skipUniqueScan and not args.unbuffered:
		logger.info( "%s unique sites (%.1f%%)" % ( unique_count, float( unique_count ) / float( site_count ) * 100.0 ) )

if __name__ == '__main__':
	run()
