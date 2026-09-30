#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

struct pt {
	double x, y;
	pt(double x, double y) : x(x), y(y) {}
	pt() {}
	bool operator<(pt p) const {
		return x < p.x - EPS || (fabs(x - p.x) <= EPS && y < p.y - EPS);
	}
};

// 7.6 Radial order of points
struct Cmp { // IMPORTANT: add const in pt operator -
	pt r;
	Cmp(pt r) : r(r) {}
	int cuad(const pt& a) const {
		if(a.x > 0 && a.y >= 0) return 0;
		if(a.x <= 0 && a.y > 0) return 1;
		if(a.x < 0 && a.y <= 0) return 2;
		if(a.x >= 0 && a.y < 0) return 3;
		assert(a.x == 0 && a.y == 0);
		return -1;
	}
	bool cmp(const pt& p1, const pt& p2) const {
		int c1 = cuad(p1), c2 = cuad(p2);
		if(c1 == c2) return p1.y * p2.x < p1.x * p2.y;
		return c1 < c2;
	}
	bool operator()(const pt& p1, const pt& p2) const {
		return cmp(pt(p1.x - r.x, p1.y - r.y),
			pt(p2.x - r.x, p2.y - r.y));
	}
};
