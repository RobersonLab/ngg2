import os

TEST_FILE = os.path.join( os.path.dirname( __file__ ), 'test.fa' )

from ngg2 import rev_comp, compile_regex_patterns, tuple_to_key, Grna, multiproc_ngg_scan

#def multiproc_ngg_scan( fasta_fname, region_tuple, regex_tuple, scan_type, max_g_content ):
# contig, start, end, strand
# def __init__( self, contig, positionStart, strand, regexMatch, inMatch=None, inStart=None )
GAGAAGACTATTTCCGTAGGAGG a plus
AAGACTATTTCCGTAGGAGGTGG b plus

def test_ngg_scan_plus_block():
	strand = '+'
	scantype = 'block'
	allow_non_g_start = False
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	assert( out == [ Grna( contig, start_pos, strand, None, 'GAGAAGACTATTTCCGTAGGAGG', 1 ) ] )
	
def test_ngg_scan_plus_exhaustive():
	strand = '+'
	scantype = 'exhaustive'
	allow_non_g_start = True
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	assert( out == [ Grna( contig, start_pos, strand, None, 'GAGAAGACTATTTCCGTAGGAGG', 0 ), Grna( contig, start_pos, strand, None, 'AAGACTATTTCCGTAGGAGGTGG', 3 ) ] )
	
def test_ngg_scan_minus_block():
	strand = '-'
	scantype = 'block'
	allow_non_g_start = False
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	assert( out == [ Grna( contig, start_pos, strand, None, 'CCTCCTACGGAAATAGTCTTCTC', 216 ) ] )

def test_ngg_scan_minus_exhaustive():
	strand = '-'
	scantype = 'exhaustive'
	allow_non_g_start = False
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	assert( out == [ Grna( contig, start_pos, strand, None, 'CCTCCTACGGAAATAGTCTTCTC', 216 ) ] )
