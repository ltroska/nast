from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "<bitset>" namespace "std" nogil:
	cdef cppclass bitset "std::bitset<7>":
		bitset() except+
		
		bitset(int) except+
				
		bool& operator[](int pos)
		
		int count()
		
		int size()
		
		bool test(int)
		
		bool any()
		
		bool none()
		
		bool all()
		
		bitset& set()
		bitset& set(int, bool)
		
		bitset& reset()
		bitset& reset(int)
		
		bitset& flip()
		bitset& flip(int)
		
		string to_string()
		
		
		
		
		
