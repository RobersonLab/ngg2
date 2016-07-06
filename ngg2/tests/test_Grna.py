from ngg2 import rev_comp

def test_rev_comp():
	assert( rev_comp( 'ACGTN' ) == 'NACGT' )
