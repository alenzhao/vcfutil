#include <cstdlib>
#include <vector>
#include <algorithm>

namespace algo {

	// partition array around pivot indexed by k
	// return the new index of the pivot
	template <typename T>
	size_t partition(T* a, size_t n, size_t k) {
		// move pivot to the first position
		std::swap(a[0], a[k]);

		// use first element as pivot
		T p = a[0];
		// i demarcates the boundary between partition "< p" and partition "> p"
		// i indexes the first element of "> p"
		size_t i = 1;
		for (size_t j = 1; j < n; ++j) {
			if (a[j] < p) {
				// swap element with the leftmost element in partition "> p"
				std::swap(a[i], a[j]);
				// increment i to include the new element in partition "< p"
				++i;
			}
		}

		// move pivot to middle, swapping with rightmost element in partition "< p"
		k = i - 1;
		std::swap(a[0], a[k]);
		return k;
	}

	// return the ith order statistic
	// post: a is rearranged
	template <typename T>
	T rselect(T* a, size_t n, size_t i) {
		if (n == 1) return a[0];

		// choose pivot p from A uniformly at random
		// and partition A around p
		size_t j = partition(a, n, std::rand() % n);

		if (j == i) return a[j];

		// recursion
		if (j > i) return rselect(a, j, i);
		if (j++ < i) return rselect(a+j, n-j, i-j);
	}

	// return the median if n is odd
	// return the greater of the middle 2 elements if n is even
	template <typename T>
	T median(T* a, size_t n) {
		return rselect(a, n, n/2);
	}

	template <typename T>
	void qsort(T* a, size_t n) {
		if (n <= 1) return;
		size_t j = partition(a, n, std::rand() % n);
		// recurse on the subarray of elements < pivot
		if (j > 0) qsort(a, j);
		// recurse on the subarray of elements > pivot
		if (++j < n) qsort(a+j, n-j);
	}

}
