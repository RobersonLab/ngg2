import os

TEST = os.path.join( os.path.dirname( __file__ ), 'test.fa' )

from ngg2 import rev_comp, compile_regex_patterns, tuple_to_key, Grna, multiproc_ngg_scan

def test_ngg_scan_plus_block():
	strand = '+'
	scantype = 'block'
	allow_non_g_start = False
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	test = [ Grna( contig, start_pos, strand, None, 'GAGAAGACTATTTCCGTAGGAGG', 0 ) ]
	
	assert( out[1] == test )
	
def test_ngg_scan_plus_exhaustive():
	strand = '+'
	scantype = 'exhaustive'
	allow_non_g_start = True
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	test = [ Grna( contig, start_pos, strand, None, 'GAGAAGACTATTTCCGTAGGAGG', 0 ), Grna( contig, start_pos, strand, None, 'AAGACTATTTCCGTAGGAGGTGG', 3 ) ]
	
	assert( out[1] == test )
	
def test_ngg_scan_minus_block():
	strand = '-'
	scantype = 'block'
	allow_non_g_start = False
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	test = [ Grna( contig, start_pos, strand, None, 'CCTCCTACGGAAATAGTCTTCTC', 217 ) ]
	
	assert( out[1] == test )

def test_ngg_scan_minus_exhaustive():
	strand = '-'
	scantype = 'exhaustive'
	allow_non_g_start = False
	contig = 'test_a'
	start_pos = 0
	
	out = multiproc_ngg_scan( TEST, ( contig, start_pos, 240, strand ), compile_regex_patterns( scantype, allow_non_g_start ), scantype, 15 )
	
	test = [ Grna( contig, start_pos, strand, None, 'CCTCCTACGGAAATAGTCTTCTC', 217 ) ]
	
	assert( out[1] == test )
