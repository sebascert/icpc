#include "template.hpp"

const int N = 1e5;
int n, a[N];
vector<ii> g[N];

lli dist[N];
void bford(int u) {
	fill(dist, dist + n, INF);
	dist[u] = 0;
	fore(i, 0, n) fore(u, 0, n) if(dist[u] != INF) for(auto [v, i]: g[u])
		dist[v] = min(dist[v], dist[u] + a[i]);
	// neg cycle all reachable from v have -INF dist
	// fore(u,0,n)if(dist[u]!=INF)for(auto[v,i]:g[u])if(dist[v]>dist[u]+a[i]);
}
