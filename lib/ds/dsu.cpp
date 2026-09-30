#include "template.hpp"

struct DSU {
	// vector<ii> st;
	vector<int> pr, rank;
	DSU(int n) : pr(n), rank(n, 1) { iota(all(pr), 0); }
	int find(int u) { return pr[u] == u ? u : pr[u] = find(pr[u]); }
	// int find(int u) { return pr[u] == u ? u : find(pr[u]); }
	void unite(int u, int v) {
		u = find(u), v = find(v);
		if(u == v) {
			// st.pb({-1, -1});
			return;
		}
		if(rank[u] < rank[v]) swap(u, v);
		pr[v] = u;
		rank[u] += rank[v];
		// st.pb({u, v});
	}
	// void rollback() {
	// 	auto [u, v] = st.back();
	// 	st.pop_back();
	// 	if(u != -1) {
	// 		pr[v] = v;
	// 		rank[u] -= rank[v];
	// 	}
	// }
};

/*
struct DSU {
	vector<int> pr, rank;
	DSU(int n) : pr(n), rank(n, 1) { iota(all(pr), 0); }
	int find(int u) { return pr[u] == u ? u : pr[u] = find(pr[u]); }
	void unite(int u, int v) {
		u = find(u), v = find(v);
		if(u == v) return;
		if(rank[u] < rank[v]) swap(u, v);
		pr[v] = u;
		rank[u] += rank[v];
	}
};
*/
