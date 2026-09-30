#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;
const double DINF = 1e18;

struct pt {
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

struct ln {
	pt p, pq;
	ln(pt p, pt q) : p(p), pq(q - p) {}
	ln() {}
	bool has(pt r) { return dist(r) <= EPS; }
	bool seghas(pt r) { return has(r) && (r - p) * (r - (p + pq)) <= EPS; }
	bool operator/(ln l) { return fabs(pq.unit() % l.pq.unit()) <= EPS; }
	bool operator==(ln l) { return *this / l && has(l.p); }
	pt operator^(ln l) {
		if(*this / l) return pt(DINF, DINF);
		pt r = l.p + l.pq * ((p - l.p) % pq / (l.pq % pq));
		return r;
	}
	double angle(ln l) { return pq.angle(l.pq); }
	int side(pt r) { return has(r) ? 0 : sgn2(pq % (r - p)); }
	pt proj(pt r) { return p + pq * ((r - p) * pq / pq.norm2()); }
	pt ref(pt r) { return proj(r) * 2 - r; }
	double dist(pt r) { return (r - proj(r)).norm(); }
	ln rot(double a) { return ln(p, p + pq.rot(a)); }
	ln rot(pt a) { return ln(p, p + pq.rot(a)); }
};

ln bisector(ln l, ln m) {
	pt p = l ^ m;
	return ln(p, p + l.pq.unit() + m.pq.unit());
}
ln bisector(pt p, pt q) { return ln((p + q) * .5, p).rot(ccw90); }

// 7.3 Circle
struct circle {
	pt o;
	double r;
	circle(pt o, double r) : o(o), r(r) {}
	circle(pt x, pt y, pt z) {
		o = bisector(x, y) ^ bisector(x, z);
		r = (o - x).norm();
	}
	bool has(pt p) { return (o - p).norm() <= r + EPS; }
	vector<pt> operator^(circle c) { // ccw
		vector<pt> s;
		double d = (o - c.o).norm();
		if(d > r + c.r + EPS || d + min(r, c.r) + EPS < max(r, c.r))
			return s;
		double x = (d * d - c.r * c.r + r * r) / (2 * d);
		double y = sqrt(r * r - x * x);
		pt v = (c.o - o) / d;
		s.push_back(o + v * x - v.rot(ccw90) * y);
		if(y > EPS) s.push_back(o + v * x + v.rot(ccw90) * y);
		return s;
	}
	vector<pt> operator^(ln l) {
		vector<pt> s;
		pt p = l.proj(o);
		double d = (p - o).norm();
		if(d - EPS > r) return s;
		if(fabs(d - r) <= EPS) {
			s.push_back(p);
			return s;
		}
		d = sqrt(r * r - d * d);
		s.push_back(p + l.pq.unit() * d);
		s.push_back(p - l.pq.unit() * d);
		return s;
	}
	vector<pt> tang(pt p) {
		double d = sqrt((p - o).norm2() - r * r);
		return *this ^ circle(p, d);
	}
	bool in(circle c) { // non strict
		double d = (o - c.o).norm();
		return d + r <= c.r + EPS;
	}
	double intertriangle(pt a, pt b) { // area of intersection with oab
		if(fabs((o - a) % (o - b)) <= EPS) return 0.0;
		vector<pt> q = {a}, w = *this ^ ln(a, b);
		if(w.size() == 2)
			for(auto p: w)
				if((a - p) * (b - p) < -EPS) q.push_back(p);
		q.push_back(b);
		if(q.size() == 4 && (q[0] - q[1]) * (q[2] - q[1]) > EPS)
			swap(q[1], q[2]);
		double s = 0;
		for(int i = 0; i < (int)q.size() - 1; ++i) {
			if(!has(q[i]) || !has(q[i + 1]))
				s += r * r * (q[i] - o).angle(q[i + 1] - o) / 2;
			else s += fabs((q[i] - o) % (q[i + 1] - o) / 2);
		}
		return s;
	}
};
