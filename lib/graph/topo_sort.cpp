#include "template.hpp"

const int N = 0;
int n;
vector<int> g[N];

// sort source to sink, for sink to source flip every edge
// if sz(sort)<n there is a cycle
vector<int> tsort() {
	queue<int> q;
	vector<int> d(2 * n, 0);
	fore(u, 0, n) for(int v: g[u]) d[v]++;
	fore(u, 0, n) if(!d[u]) q.push(u);
	vector<int> r;
	while(!q.empty()) {
		int u = q.front();
		q.pop();
		r.pb(u);
		for(int v: g[u]) {
			d[v]--;
			if(!d[v]) q.push(v);
		}
	}
	return r;
}
