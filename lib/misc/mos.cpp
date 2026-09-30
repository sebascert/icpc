#include "template.hpp"

const int Q = 0;

void init();
void add(int i);
void rem(int i);
lli ans();

struct qu {
	int l, r, id;
};
int n, sq, nq; // array size, sqrt(array size), #queries
qu qs[Q];
lli qq[Q];
bool qcomp(const qu& a, const qu& b) {
	if(a.l / sq != b.l / sq) return a.l < b.l;
	return (a.l / sq) & 1 ? a.r < b.r : a.r > b.r;
}
void mos() {
	sq = sqrt(n) + .5;
	sort(qs, qs + nq, qcomp);
	init();
	int l = 0, r = 0;
	fore(i, 0, nq) {
		qu q = qs[i];
		while(l > q.l) add(--l);
		while(r < q.r) add(r++);
		while(l < q.l) rem(l++);
		while(r > q.r) rem(--r);
		qq[q.id] = ans();
	}
}
