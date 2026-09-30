#include "template.hpp"
#include "ds/dsu.hpp"

const int N = 0;
int n, m;
vector<ii> g[N];

lli mst() {
	vector<vi> e;
	fore(u, 0, n) for(auto [v, w]: g[u]) e.pb({w, u, v});
	sort(all(e));
	DSU dsu(n);
	lli x = 0;
	fore(i, 0, m) {
		lli u = e[i][1], v = e[i][2], w = e[i][0];
		if(dsu.unite(u, v)) x += w;
	}
	return x;
}
