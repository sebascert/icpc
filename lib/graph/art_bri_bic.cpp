#include "template.hpp"
#include "ds/dsu.hpp"

const int N = 0 + 1;
int n, m;
vector<int> t[N];

DSU dsu(0);
int d[N], art[N];
int dfs(int u, int p) {
	int dp = d[u] = d[p] + 1;
	for(int v: t[u]) {
		if(v == p) continue;
		if(!d[v]) { // span edge
			dp = min(dp, dfs(v, u));
		} else if(d[v] < d[u]) { // back edge
			dp = min(dp, d[v]);
		}
	}
	art[p] += dp >= d[p];
	if(p != n) {
		if(dp <= d[u]) { // bridge (u,p)
		} else {
			dsu.unite(u, p);
		}
	} else art[u]--; // root
	return dp;
}
void art_bri_bic() {
	d[n] = 0;
	dsu = DSU(n);
	dfs(0, n);
}
