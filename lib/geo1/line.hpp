#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;
const double DINF = 1e18;

struct pt { // for 3D add z coordinate
	double x, y;
	pt(double x, double y) : x(x), y(y) {}
	pt() {}
	double norm2() { return *this * *this; }
	double norm() { return sqrt(norm2()); }
	bool operator==(pt p) {
		return fabs(x - p.x) <= EPS && fabs(y - p.y) <= EPS;
	}
	pt operator+(pt p) { return pt(x + p.x, y + p.y); }
	pt operator-(pt p) { return pt(x - p.x, y - p.y); }
	pt operator*(double t) { return pt(x * t, y * t); }
	pt operator/(double t) { return pt(x / t, y / t); }
	double operator*(pt p) { return x * p.x + y * p.y; }
	double angle(pt p) { return acos((*this * p) / (norm() * p.norm())); }
	pt unit() { return *this / norm(); }
	double operator%(pt p) { return x * p.y - y * p.x; }
	bool operator<(pt p) const {
		return x < p.x - EPS || (fabs(x - p.x) <= EPS && y < p.y - EPS);
	}
	bool left(pt p, pt q) { return (q - p) % (*this - p) > EPS; }
	pt rot(pt r) { return pt((*this) % r, (*this) * r); }
	pt rot(double a) { return rot(pt(sin(a), cos(a))); }
};

pt ccw90(1, 0);
pt cw90(-1, 0);

static inline int sgn2(double x) { return (x < -EPS) ? -1 : (x > EPS ? 1 : 0); }

// 7.2 Line
// Needs EPS, INF and DINF
struct ln {
	pt p, pq;
	ln(pt p, pt q) : p(p), pq(q - p) {}
	ln() {}
	bool has(pt r) { return dist(r) <= EPS; }
	bool seghas(pt r) { return has(r) && (r - p) * (r - (p + pq)) <= EPS; }
	// bool operator /(ln l){return (pq.unit()^l.pq.unit()).norm()<=EPS;}
	// 3D
	bool operator/(ln l) {
		return fabs(pq.unit() % l.pq.unit()) <= EPS;
	} // 2D
	bool operator==(ln l) { return *this / l && has(l.p); }
	pt operator^(ln l) { // intersection
		if(*this / l) return pt(DINF, DINF);
		pt r = l.p + l.pq * ((p - l.p) % pq / (l.pq % pq));
		return r;
	}
	double angle(ln l) { return pq.angle(l.pq); }
	int side(pt r) { return has(r) ? 0 : sgn2(pq % (r - p)); } // 2D
	pt proj(pt r) { return p + pq * ((r - p) * pq / pq.norm2()); }
	pt ref(pt r) { return proj(r) * 2 - r; }
	double dist(pt r) { return (r - proj(r)).norm(); }
	// double dist(ln l){ // only 3D
	// if(*this/l)return dist(l.p);
	// return abs((l.p-p)*(pq^l.pq))/(pq^l.pq).norm();
	// }
	ln rot(double a) { return ln(p, p + pq.rot(a)); } // 2D
	ln rot(pt a) { return ln(p, p + pq.rot(a)); }
};

ln bisector(ln l, ln m) { // angle bisector
	pt p = l ^ m;
	return ln(p, p + l.pq.unit() + m.pq.unit());
}
ln bisector(pt p, pt q) { // segment bisector (2D)
	return ln((p + q) * .5, p).rot(ccw90);
}
