from ngg2 import tuple_to_key

def test_tuple_to_key():
	assert( tuple_to_key( 1, 'brown', '1', 15, 'cloud' ) == '1_brown_1_15_cloud' )
