#include "template.hpp"

const int N = 0;
lli n, m, a[N];
vector<ii> g[N];

lli dist[N];
void dijkstra(int u) {
	memset(dist, -1, sizeof(dist));
	priority_queue<ii> q;
	dist[u] = 0;
	q.push({u, 0});
	while(!q.empty()) {
		auto [u, d] = q.top();
		q.pop();
		if(dist[u] != d) continue;
		for(auto [v, i]: g[u]) {
			lli dv = d + a[i];
			if(dist[v] < 0 || dist[v] > dv)
				dist[v] = dv, q.push({v, dv});
		}
	}
}
