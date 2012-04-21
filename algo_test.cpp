#include <iostream>
using namespace std;

#include "algo.hpp"

void print_array(int a[], size_t n) {
	for (size_t i = 0; i < n; ++i) {
		cout << a[i] << ' ';
	}
	cout << endl;
}

int main() {

	const size_t n = 5, r = 4;
	int a[n] = {9, 3, 2, 2, -4};

	print_array(a, n);
	
	cout << "rselect[" << r << "] : " << algo::rselect(a, n, r-1) << endl;

	cout << "before qsort:" << endl;

	print_array(a, n);

	cout << "after qsort:" << endl;

	algo::qsort(a, n);

	print_array(a, n);

	return 0;
}

