#pragma once
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

struct pt {
	double x, y;
	pt(double x, double y) : x(x), y(y) {}
	pt() {}
	pt operator-(pt p) const { return pt(x - p.x, y - p.y); }
};

// 7.6 Radial order of points (used here)
struct Cmp {
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

// 7.8 Dual from planar graph
const int MAXN = 1;
vector<int> g[MAXN];
int n; // input graph (must be connected)
vector<int> gd[MAXN];
int nd;			 // output graph
vector<int> nodes[MAXN]; // nodes delimiting region (in CW order)
map<pair<int, int>, int> ps, es;

void get_dual(vector<pt> p) { // p: points corresponding to nodes
	ps.clear();
	es.clear();
	for(int x = 0; x < n; ++x) {
		Cmp pc(p[x]);
		auto comp = [&](int a, int b) { return pc(p[a], p[b]); };
		sort(g[x].begin(), g[x].end(), comp);
		for(int i = 0; i < (int)g[x].size(); ++i) ps[{x, g[x][i]}] = i;
	}
	nd = 0;
	for(int xx = 0; xx < n; ++xx)
		for(auto yy: g[xx])
			if(!es.count({xx, yy})) {
				int x = xx, y = yy;
				gd[nd].clear();
				nodes[nd].clear();
				while(!es.count({x, y})) {
					es[{x, y}] = nd;
					nodes[nd].push_back(y);
					int z = g[y][(ps[{y, x}] + 1) %
						     g[y].size()];
					x = y;
					y = z;
				}
				nd++;
			}
	for(auto p: es) {
		pair<int, int> q = {p.first.second, p.first.first};
		assert(es.count(q));
		if(es[q] != p.second) gd[p.second].push_back(es[q]);
	}
	for(int i = 0; i < nd; ++i) {
		sort(gd[i].begin(), gd[i].end());
		gd[i].erase(unique(gd[i].begin(), gd[i].end()), gd[i].end());
	}
}
