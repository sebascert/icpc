#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;
const double INF = 1e18;

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

static inline int sgn(double x) { return x < -EPS ? -1 : x > EPS; }
static inline int sgn2(double x) { return x < -EPS ? -1 : x > EPS ? 1 : 0; }

struct ln {
	pt p, pq;
	ln(pt p, pt q) : p(p), pq(q - p) {}
	ln() {}
	bool has(pt r) { return dist(r) <= EPS; }
	bool seghas(pt r) { return has(r) && (r - p) * (r - (p + pq)) <= EPS; }
	bool operator/(ln l) { return fabs(pq.unit() % l.pq.unit()) <= EPS; }
	bool operator==(ln l) { return *this / l && has(l.p); }
	pt operator^(ln l) {
		if(*this / l) return pt(INF, INF);
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

struct circle {
	pt o;
	double r;
	circle(pt o, double r) : o(o), r(r) {}
	circle(pt x, pt y, pt z) {
		o = ln(x, y) ^ ln(x, z);
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

struct pol {
	int n;
	vector<pt> p;
	pol() {}
	pol(vector<pt> _p) {
		p = _p;
		n = (int)p.size();
	}
	double area() {
		double r = 0.0;
		for(int i = 0; i < n; ++i) r += p[i] % p[(i + 1) % n];
		return fabs(r) / 2; // negative if CW, positive if CCW
	}
	pt centroid() { // barycenter
		pt r(0, 0);
		double t = 0;
		for(int i = 0; i < n;) {
			r = r +
			    (p[i] + p[(i + 1) % n]) * (p[i] % p[(i + 1) % n]);
			t += p[i] % p[(i + 1) % n];
			++i;
		}
		return r / t / 3;
	}
	bool has(pt q) { // O(n)
		for(int i = 0; i < n; ++i)
			if(ln(p[i], p[(i + 1) % n]).seghas(q)) return true;
		int cnt = 0;
		for(int i = 0; i < n; ++i) {
			int j = (i + 1) % n;
			int k = sgn((q - p[j]) % (p[i] - p[j]));
			int u = sgn(p[i].y - q.y), v = sgn(p[j].y - q.y);
			if(k > 0 && u < 0 && v >= 0) cnt++;
			if(k < 0 && v < 0 && u >= 0) cnt--;
		}
		return cnt != 0;
	}
	void normalize() { // call before haslog, remove collinear first
		if(p[2].left(p[0], p[1])) reverse(p.begin(), p.end());
		int pi = min_element(p.begin(), p.end()) - p.begin();
		vector<pt> s(n);
		for(int i = 0; i < n; ++i) s[i] = p[(pi + i) % n];
		p.swap(s);
	}
	bool haslog(pt q) { // O(log(n)) only CONVEX. Call normalize first
		if(q.left(p[0], p[1]) || q.left(p.back(), p[0])) return false;
		int a = 1, b = (int)p.size() - 1;
		while(b - a > 1) {
			int c = (a + b) / 2;
			if(!q.left(p[0], p[c])) a = c;
			else b = c;
		}
		return !q.left(p[a], p[a + 1]);
	}
	pt farthest(pt v) { // O(log(n)) only CONVEX
		if(n < 10) {
			int k = 0;
			for(int i = 1; i < n; ++i)
				if(v * (p[i] - p[k]) > EPS) k = i;
			return p[k];
		}
		if(n == (int)p.size()) p.push_back(p[0]);
		pt a = p[1] - p[0];
		int s = 0, e = n, ua = v * a > EPS;
		if(!ua && v * (p[n - 1] - p[0]) <= EPS) return p[0];
		while(1) {
			int m = (s + e) / 2;
			pt c = p[m + 1] - p[m];
			int uc = v * c > EPS;
			if(!uc && v * (p[m - 1] - p[m]) <= EPS) return p[m];
			if(ua && (!uc || v * (p[s] - p[m]) > EPS)) e = m;
			else if(ua || uc || v * (p[s] - p[m]) >= -EPS)
				s = m, a = c, ua = uc;
			else e = m;
			assert(e > s + 1);
		}
	}
	double intercircle(circle c) {
		double r = 0;
		for(int i = 0; i < n; ++i) {
			int j = (i + 1) % n;
			double w = c.intertriangle(p[i], p[j]);
			if((p[j] - c.o) % (p[i] - c.o) > 0) r += w;
			else r -= w;
		}
		return fabs(r);
	}
	double callipers() {  // square distance of most distant points
		double r = 0; // prereq: convex, ccw, NO COLLINEAR POINTS
		for(int i = 0, j = n < 2 ? 0 : 1; i < j; ++i) {
			for(;; j = (j + 1) % n) {
				r = max(r, (p[i] - p[j]).norm2());
				if((p[(i + 1) % n] - p[i]) %
						(p[(j + 1) % n] - p[j]) <=
					EPS)
					break;
			}
		}
		return r;
	}
};

// Dynamic convex hull trick
vector<pol> w;
void add(pt q) { // add(q), O(log^2(n))
	vector<pt> p = {q};
	while(!w.empty() && (int)(w.back().p).size() < 2 * (int)p.size()) {
		for(pt v: w.back().p) p.push_back(v);
		w.pop_back();
	}
	w.push_back(pol(chull(p)));
}
long long query(pt v) { // max(q*v:q in w), O(log^2(n))
	long long r = -INF;
	for(auto& p: w) r = max(r, (long long)(p.farthest(v) * v));
	return r;
}
