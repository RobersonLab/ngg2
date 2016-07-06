import re
import regex

from ngg2 import compile_regex_patterns

def test_compile_regex_patterns_block():
	compile_block_plus_noncanon  = re.compile( "[ACGT]{18}GG[ACGT]{1}GG" )
	compile_block_minus_noncanon = re.compile( "CC[ACGT]{1}CC[ACGT]{18}" )
	compile_block_plus_canon = re.compile( "G[ACGT]{17}GG[ACGT]{1}GG" )
	compile_block_minus_canon = re.compile( "CC[ACGT]{1}CC[ACGT]{17}C" )
	
	assert( compile_regex_patterns( 'block', False ) == ( compile_block_plus_canon, compile_block_minus_canon ) )
	assert( compile_regex_patterns( 'block', True ) == ( compile_block_plus_noncanon, compile_block_minus_noncanon ) )
	
def test_compile_regex_patterns_exhaustive():
	compile_exhaustive_plus_noncanon  = regex.compile( "[ACGT]{18}GG[ACGT]{1}GG" )
	compile_exhaustive_minus_noncanon = regex.compile( "CC[ACGT]{1}CC[ACGT]{18}" )
	compile_exhaustive_plus_canon = regex.compile( "G[ACGT]{17}GG[ACGT]{1}GG" )
	compile_exhaustive_minus_canon = regex.compile( "CC[ACGT]{1}CC[ACGT]{17}C" )
	
	assert( compile_regex_patterns( 'exhaustive', False ) == ( compile_exhaustive_plus_canon, compile_exhaustive_minus_canon ) )
	assert( compile_regex_patterns( 'exhaustive', True ) == ( compile_exhaustive_plus_noncanon, compile_exhaustive_minus_noncanon ) )
