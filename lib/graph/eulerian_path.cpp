#include "template.hpp"

const int N = 0;

// uncomment commented code for undirected
struct edge {
	int y;
	// list<edge>::iterator rev;
	edge(int y) : y(y) {}
};
list<edge> g[N];
void add_edge(int a, int b) {
	g[a].push_front(edge(b));
	// auto ia = g[a].begin();
	// g[b].push_front(edge(a));
	// auto ib = g[b].begin();
	// ia->rev = ib;
	// ib->rev = ia;
}
vector<int> p;
void go(int x) {
	while(g[x].size()) {
		int y = g[x].front().y;
		// g[y].erase(g[x].front().rev);
		g[x].pop_front();
		go(y);
	}
	p.push_back(x);
}
vector<int> get_path(int x) {
	p.clear();
	go(x);
	reverse(p.begin(), p.end());
	return p;
}
