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
	pt operator-(pt p) const { return pt(x - p.x, y - p.y); }
	double operator%(pt p) const { return x * p.y - y * p.x; }
	bool left(pt p, pt q) { return (q - p) % (*this - p) > EPS; }
};

// 7.7 Convex hull
vector<pt> chull(vector<pt> p) {
	if(p.size() < 3) return p;
	vector<pt> r;
	sort(p.begin(), p.end());		 // first x, then y
	for(int i = 0; i < (int)p.size(); ++i) { // lower hull
		while(r.size() >= 2 && r.back().left(r[r.size() - 2], p[i]))
			r.pop_back();
		r.push_back(p[i]);
	}
	r.pop_back();
	int k = (int)r.size();
	for(int i = (int)p.size() - 1; i >= 0; --i) { // upper hull
		while((int)r.size() >= k + 2 &&
			r.back().left(r[(int)r.size() - 2], p[i]))
			r.pop_back();
		r.push_back(p[i]);
	}
	r.pop_back();
	return r;
}
