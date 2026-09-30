#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;
const double INF = 1e18;

struct pt {
	double x, y;
	pt(double x, double y) : x(x), y(y) {}
	pt() {}
	pt operator+(pt p) const { return pt(x + p.x, y + p.y); }
	pt operator-(pt p) const { return pt(x - p.x, y - p.y); }
	pt operator*(double t) const { return pt(x * t, y * t); }
	pt operator/(double t) const { return pt(x / t, y / t); }
	double operator*(pt p) const { return x * p.x + y * p.y; }
	double operator%(pt p) const { return x * p.y - y * p.x; }
	double norm2() const { return (*this) * (*this); }
	double norm() const { return sqrt(norm2()); }
	pt unit() const { return *this / norm(); }
	pt rot(pt r) const { return pt((*this) % r, (*this) * r); }
	pt rot(double a) const { return rot(pt(sin(a), cos(a))); }
	bool left(pt p, pt q) const { return (q - p) % (*this - p) > EPS; }
	bool operator<(pt p) const {
		return x < p.x - EPS || (fabs(x - p.x) <= EPS && y < p.y - EPS);
	}
};

pt ccw90(1, 0);

struct ln {
	pt p, pq;
	ln(pt p, pt q) : p(p), pq(q - p) {}
	ln() {}
	pt operator^(ln l) const {
		double t = (l.p - p) % l.pq / (pq % l.pq);
		return p + pq * t;
	}
};

// 7.9 Halfplane intersection
struct halfplane : public ln {
	double angle;
	halfplane() {}
	halfplane(pt a, pt b) : p(a), pq(b - a), angle(atan2(pq.y, pq.x)) {}
	bool operator<(halfplane const& l) const { return angle < l.angle; }
	bool out(pt q) { return pq % (q - p) > EPS; }
};

vector<pt> intersect(vector<halfplane> h) {
	sort(h.begin(), h.end());
	deque<halfplane> q;
	deque<pt> p;
	for(int i = 0; i < (int)h.size(); ++i) {
		while(q.size() > 1 && h[i].out(p.back()))
			q.pop_back(), p.pop_back();
		while(q.size() > 1 && h[i].out(p.front()))
			q.pop_front(), p.pop_front();
		if(!q.empty() && fabs((q.back().pq) % (h[i].pq)) <= EPS) {
			if((h[i].pq % (q.back().p - h[i].p)) > EPS)
				q.pop_back();
			else continue;
		}
		if(!q.empty()) p.push_back(q.back() ^ h[i]);
		q.push_back(h[i]);
	}
	while(q.size() > 2 && q.front().out(p.back()))
		q.pop_back(), p.pop_back();
	while(q.size() > 2 && q.back().out(p.front()))
		q.pop_front(), p.pop_front();
	if(q.size() < 3) return {};
	p.push_back(q.back() ^ q.front());
	return vector<pt>(p.begin(), p.end());
}

vector<pt> intersectRandomized(vector<halfplane> h) {
	random_shuffle(h.begin(), h.end());
	deque<halfplane> dq;
	deque<pt> poly;
	for(int i = 0; i < (int)h.size(); ++i) {
		while(dq.size() >= 2 && h[i].out(poly.back()))
			dq.pop_back(), poly.pop_back();
		while(dq.size() >= 2 && h[i].out(poly.front()))
			dq.pop_front(), poly.pop_front();
		if(!dq.empty() && fabs(dq.back().pq % h[i].pq) <= EPS) {
			if(h[i].pq % (dq.back().p - h[i].p) > EPS)
				dq.pop_back();
			else continue;
		}
		if(!dq.empty()) poly.push_back(dq.back() ^ h[i]);
		dq.push_back(h[i]);
	}
	while(dq.size() >= 3 && dq.front().out(poly.back()))
		dq.pop_back(), poly.pop_back();
	while(dq.size() >= 3 && dq.back().out(poly.front()))
		dq.pop_front(), poly.pop_front();
	if(dq.size() < 3) return {};
	poly.push_back(dq.back() ^ dq.front());
	return vector<pt>(poly.begin(), poly.end());
}
