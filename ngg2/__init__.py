#!/usr/bin/env python

"""Tool to find 3'GG CRISPR/Cas9 sites for high-efficiency genome-editing"""

##########
# import #
##########
import pyfaidx
import regex
import re
import sys
import time

####################
# Version and name #
####################
__script_path__ = sys.argv[0]
__script_name__ = __script_path__.split( '/' )[-1].split( '\\' )[-1]
__version__ = 'v1.3.0'

########
# fxns #
########
def rev_comp( seq ):
	"""DNA reverse complement
	"""
	
	nuc_dict = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N' }
	return ''.join( [ nuc_dict[c] for c in seq.upper()[::-1] ] )

def compile_regex_patterns( scan_type, allow_noncanonical ):
	"""Compile appropriate regular expressions to find + or - strand 3'GG sites
	"""
	
	if type( allow_noncanonical ) is not bool:
		raise( TypeError( "In compile_regex_patterns, allow_noncanonical must be a bool!" ) )
	elif scan_type not in [ 'exhaustive', 'block' ]:
		raise( ValueError( "In compile_regex_patterns, only exhaustive and block are valid options!" ) )
	
	if scan_type == 'exhaustive':
		if allow_noncanonical == False:
			sense_regex = regex.compile( 'G[ACGT]{17}GG[ACGT]{1}GG' )
			antisense_regex = regex.compile( 'CC[ACGT]{1}CC[ACGT]{17}C' )
		elif allow_noncanonical == True:
			sense_regex = regex.compile( '[ACGT]{18}GG[ACGT]{1}GG' )
			antisense_regex = regex.compile( 'CC[ACGT]{1}CC[ACGT]{18}' )
	
	elif scan_type == 'block':
		if allow_noncanonical == False:
			sense_regex = re.compile( 'G[ACGT]{17}GG[ACGT]{1}GG' )
			antisense_regex = re.compile( 'CC[ACGT]{1}CC[ACGT]{17}C' )
		elif allow_noncanonical == True:
			sense_regex = re.compile( '[ACGT]{18}GG[ACGT]{1}GG' )
			antisense_regex = re.compile( 'CC[ACGT]{1}CC[ACGT]{18}' )
	
	return ( sense_regex, antisense_regex )

def tuple_to_key( tuple ):
	return '_'.join( [ str(x) for x in tuple ] )

###########
# classes #
###########
class Grna:
	"""
	contig = name of the contig from the FASTA file processed
	positionStart = start position of region scanned*
	strand = strand where the match was found with regard to FASTA file
	regexMatch = the regex match object from the compiled search
	
	The positionStart is expected to be 0-base. Very important!
	"""
	
	def __init__( self, contig, positionStart, strand, regexMatch, inMatch=None, inStart=None ):
		if strand not in [ '+', '-' ]:
			raise( ValueError( "In Grna object, only '+' and '-' are acceptable strand designations" ) )
		
		if regexMatch is not None:		
			if strand == '+':
				self.seq = regexMatch.group()[:20]
				self.pam = regexMatch.group()[20:]
				self.start = positionStart + regexMatch.start() + 1 # switch to 1-base
				self.end = self.start + len( self.seq ) - 1
			elif strand == '-':
				self.seq = rev_comp( regexMatch.group()[3:] )
				self.pam = rev_comp( regexMatch.group()[:3] )
				self.end = positionStart + regexMatch.start() + 3 + 1 # switch to 1-base
				self.start = self.end + len( self.seq ) - 1
		else:
			if strand == '+':
				self.seq = inMatch[:20]
				self.pam = inMatch[20:]
				self.start = positionStart + inStart + 1
				self.end = self.start + len( self.seq ) - 1
			elif strand == '-':
				self.seq = rev_comp( inMatch[3:] )
				self.pam = rev_comp( inMatch[:3] )
				self.end = positionStart + inStart + 3 + 1
				self.start = self.end + len( self.seq ) - 1
		
		self.gCount = self.seq.count( 'G' )
		self.contig = contig
		self.strand = strand
		self.gStart = True if self.seq[0] == 'G' else False
	
	def __str__( self ):
		return "%s,%s,%s,%s,%s,%s,%s" % ( self.contig, self.start, self.end, self.seq, self.pam, self.strand, self.gStart )
		   
	def __repr__( self ):
		return str( self )
		
	def __eq__( self, other ):
		if self.seq == other.seq and self.pam == other.pam and self.start == other.start and self.end == other.end and self.strand == other.strand and self.contig == other.contig:
			return True
		else:
			return False

######################
# fxns using classes #
######################
def unbuffered_ngg_scan( fasta_fname, output_fname, region_list, regex_tuple, scan_type, max_g_content ):
	"""
	fasta_fname = path to the FASTA file to process
	output_fname = name of the output file you want to generate
	region_list = tuple with region info for scan (contig, start*, end, strand)
	regex_tuple = the + and - strand compiled regular expressions to search for
	scan_type = block or exhaustive
	max_g_content = max G fraction of bases in motif to keep
	
	*region_list is expected in 0-start begin and 1-start end from pyfaidx. If you generate the tuple yourself, make sure it is in that format.
	"""
	
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
	"""
	fasta_fname = path to the FASTA file to process
	region_tuple = tuple with region info for scan (contig, start*, end, strand)
	regex_tuple = the + and - strand compiled regular expressions to search for
	scan_type = block or exhaustive
	max_g_content = max G fraction of bases in motif to keep
	
	*region_tuple is expected in 0-start begin and 1-start end from pyfaidx. If you generate the tuple yourself, make sure it is in that format.
	"""
	
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
	pass
	