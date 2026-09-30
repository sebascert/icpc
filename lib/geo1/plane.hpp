#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;
const double DINF = 1e18;

// 7.5 Plane (3D)
// This section assumes a 3D point with cross product.
struct pt {
	double x, y, z;
	pt(double x, double y, double z) : x(x), y(y), z(z) {}
	pt() : x(0), y(0), z(0) {}
	pt operator+(pt p) { return pt(x + p.x, y + p.y, z + p.z); }
	pt operator-(pt p) { return pt(x - p.x, y - p.y, z - p.z); }
	pt operator*(double t) { return pt(x * t, y * t, z * t); }
	pt operator/(double t) { return pt(x / t, y / t, z / t); }
	double operator*(pt p) { return x * p.x + y * p.y + z * p.z; }
	pt operator^(pt p) {
		return pt(y * p.z - z * p.y, z * p.x - x * p.z,
			x * p.y - y * p.x);
	}
	double norm2() { return *this * *this; }
	double norm() { return sqrt(norm2()); }
	pt unit() { return *this / norm(); }
};

struct ln {
	pt p, pq;
	ln(pt p, pt q) : p(p), pq(q - p) {}
	ln() {}
};

struct plane {
	pt a, n; // n: normal unit vector
	plane(pt a, pt b, pt c) : a(a), n(((b - a) ^ (c - a)).unit()) {}
	plane() {}
	bool has(pt p) { return fabs((p - a) * n) <= EPS; }
	double angle(plane w) { return acos(n * w.n); }
	double dist(pt p) { return fabs((p - a) * n); }
	pt proj(pt p) {
		inter(ln(p, p + n), p);
		return p;
	}
	bool inter(ln l, pt& r) {
		double x = n * (l.p + l.pq - a), y = n * (l.p - a);
		if(fabs(x - y) <= EPS) return false;
		r = (l.p * x - (l.p + l.pq) * y) / (x - y);
		return true;
	}
	bool inter(plane w, ln& r) {
		pt nn = n ^ w.n;
		pt v = n ^ nn;
		double d = w.n * v;
		if(fabs(d) <= EPS) return false;
		pt p = a + v * (w.n * (w.a - a) / d);
		r = ln(p, p + nn);
		return true;
	}
};
