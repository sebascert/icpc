#include "template.hpp"

const int N = 0 + 1, sqN = 0;
lli sq[sqN];
int a[N], who[N];
int block;
int n;

template <class T> T query(int l, int r) {
	T x = {};
	while(l <= r) {
		if(l % block == 0 && l + block - 1 <= r) {
			x += sq[who[l]];
		} else {
			x += a[l++];
		}
	}
	return x;
}
void do_sqrt_decomposition() {
	block = sqrt(n);
	fore(i, 0, n) {
		who[i] = i / block;
		sq[who[i]] += a[i];
	}
}
