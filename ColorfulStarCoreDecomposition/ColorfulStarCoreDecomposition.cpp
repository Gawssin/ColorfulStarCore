#include <stdlib.h>
#include <cstdio>
#include <stdbool.h>
#include <string.h>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <cstdarg>
#include "heapLLU.h"

#define NLINKS 1000000
using namespace std;

int* loc, * tmploc, cvCnt = 0;
long long kCCN = 0, *tmpCnt, * cntSub;


int* uadj, * uadjt, * Merge, * mark, * oloc;

keyvalueLLU* cv;
class edge
{
public:
	int s;
	int t;
};

class iddeg
{
public:
	int id;
	int degree;
};

class Graph
{
public:
	Graph();
	Graph(const Graph& obj);
	~Graph();
	void readedgelist(string edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	Graph* mksub(int, ...);
	//Graph* mksubMark(int*, int, int*);
	int color(int*);
	void kClique(int, long long*, long long*);
	void kCliqueCount(int, long long*, int*, int*, int*, int*, int*, int**, int**, long long*);
	//bool isEdge(int, int);

	int n;
	int e;
	int maxDeg;
	edge* edges;

	int* deg;
	int* cd;
	int* adj;
	int* coreRank;			//increasing core number order
	int* coreNum;			//coreNum[i] is the core number of node i.
	int* bin;
};

inline int max3(int a, int b, int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

Graph::Graph(void) {}
Graph::~Graph(void)
{
	if (deg != NULL) delete[] deg;
	if (cd != NULL) delete[] cd;
	if (adj != NULL) delete[] adj;
	if (coreRank != NULL) delete[] coreRank;
	if (coreNum != NULL) delete[] coreNum;
	if (bin != NULL) delete[] bin;
}
Graph::Graph(const Graph& obj)
{
	n = obj.n, e = obj.e, maxDeg = obj.maxDeg, edges = obj.edges;
	edges = obj.edges;

	if (deg != NULL) delete[] deg;
	if (obj.deg != NULL) deg = new int[n], memcpy(deg, obj.deg, n * sizeof(int));

	if (cd != NULL) delete[] cd;
	if (obj.cd != NULL) cd = new int[n + 1], memcpy(cd, obj.cd, (n + 1) * sizeof(int));

	if (adj != NULL) delete[] adj;
	if (obj.adj != NULL) adj = new int[2 * e], memcpy(adj, obj.adj, 2 * e * sizeof(int));

	if (coreRank != NULL) delete[] coreRank;
	if (obj.coreRank != NULL) coreRank = new int[n], memcpy(coreRank, obj.coreRank, n * sizeof(int));

	if (coreNum != NULL) delete[] coreNum;
	if (obj.coreNum != NULL) coreNum = new int[n], memcpy(coreNum, obj.coreNum, n * sizeof(int));

	if (bin != NULL) delete[] bin;
	if (obj.bin != NULL) bin = new int[maxDeg + 2], memcpy(bin, obj.bin, (maxDeg + 2) * sizeof(int));
}

void Graph::readedgelist(string edgelist)
{
	ifstream file;
	file.open(edgelist);
	file >> n >> e;
	edges = new edge[e];
	e = 0;
	while (file >> edges[e].s >> edges[e].t) e++;
	file.close();
}

void Graph::mkGraph()
{
	deg = new int[n]();
	cd = new int[n + 1];
	adj = new int[2 * e];
	maxDeg = 0;
	for (int i = 0; i < e; i++)
	{
		deg[edges[i].s]++;
		deg[edges[i].t]++;
		maxDeg = max3(maxDeg, deg[edges[i].s], deg[edges[i].t]);
	}
	cd[0] = 0;
	for (int i = 1; i < n + 1; i++)
	{
		cd[i] = cd[i - 1] + deg[i - 1];
		deg[i - 1] = 0;
	}

	for (int i = 0; i < e; i++)
	{
		adj[cd[edges[i].s] + deg[edges[i].s]++] = edges[i].t;
		adj[cd[edges[i].t] + deg[edges[i].t]++] = edges[i].s;
	}

	for (int i = 0; i < n; i++) sort(adj + cd[i], adj + cd[i] + deg[i]);
}

bool cmp(const pair<int, int>& a, const pair<int, int>& b)
{
	return a.second > b.second;
}
bool IGCmp(const iddeg& a, const iddeg& b)
{
	return a.degree == b.degree ? (a.id < b.id) : (a.degree > b.degree);
}

bool Graph::isEdge(int a, int b)
{
	if (deg[a] > deg[b]) a = a ^ b, b = a ^ b, a = a ^ b;
	for (int i = cd[a]; i < cd[a] + deg[a]; i++)
		if (adj[i] == b) return true;
	return false;
}
int Graph::outLargeClique()
{
	int CSize = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		int id = coreRank[i];
		if (coreNum[id] >= CSize)
		{
			pair<int, int>* SCore = new pair<int, int>[deg[id]];
			//int *S = new int[deg[id]], cnt = 0;
			int cnt = 0, ind = 0;
			for (int j = cd[id]; j < cd[id] + deg[id]; j++)
				if (coreNum[adj[j]] >= CSize)
				{
					SCore[cnt].first = adj[j];
					SCore[cnt].second = coreNum[adj[j]];
					cnt++;
				}
			//S[cnt++] = adj[j];
			sort(SCore, SCore + cnt, cmp);
			//sort(S, S + cnt, cmp);
			int* C = new int[deg[id]];
			for (int j = 0; j < cnt; j++)
			{
				int flag = 1;
				for (int k = 0; k < ind; k++)
				{
					if (isEdge(SCore[j].first, C[k]) == false)
					{
						flag = 0;
						break;
					}
				}
				if (flag) C[ind++] = SCore[j].first;
			}
			ind++;	//node "id" ?
			if (ind > CSize) CSize = ind;
		}
	}
	return CSize;
}

void Graph::coreDecomposition()
{
	bin = new int[maxDeg + 2]();

	for (int i = 0; i < n; i++)
		bin[deg[i]]++;

	int lastBin = bin[0], nowBin;
	bin[0] = 0;
	for (int i = 1; i <= maxDeg; i++)
	{
		nowBin = lastBin + bin[i - 1];
		lastBin = bin[i];
		bin[i] = nowBin;
	}
	int* vert = new int[n](), * pos = new int[n](), * tmpDeg = new int[n]();
	for (int i = 0; i < n; i++)
	{
		pos[i] = bin[deg[i]];

		vert[bin[deg[i]]++] = i;
		tmpDeg[i] = deg[i];
	}

	bin[0] = 0;
	for (int i = maxDeg; i >= 1; i--)
	{
		bin[i] = bin[i - 1];
	}

	//int core = 0;
	int* cNum = new int[n];
	//int *cNum = (int *)malloc(g->n * sizeof(int));
	for (int i = 0; i < n; i++)
	{
		int id = vert[i], nbr, binFrontId;
		//if (i == bin[core + 1]) ++core;
		cNum[id] = tmpDeg[id];
		for (int i = cd[id]; i < cd[id] + deg[id]; i++)
		{
			nbr = adj[i];

			if (tmpDeg[nbr] > tmpDeg[id])
			{
				binFrontId = vert[bin[tmpDeg[nbr]]];
				if (binFrontId != nbr)
				{

					pos[binFrontId] = pos[nbr];
					pos[nbr] = bin[tmpDeg[nbr]];
					vert[bin[tmpDeg[nbr]]] = nbr;
					vert[pos[binFrontId]] = binFrontId;

				}
				bin[tmpDeg[nbr]]++;
				tmpDeg[nbr]--;

			}

		}

	}

	coreNum = cNum;

	coreRank = vert;

	delete[] tmpDeg;
	delete[] pos;
}

Graph* Graph::mksub(int argCnt, ...)//(int *nodes, int NodeNum)//
{
	Graph* sg = new Graph;
	int* Mark = NULL;
	va_list ap;
	va_start(ap, argCnt);
	//int *sg = va_arg(ap,int *);

	int* nodes = va_arg(ap, int*);
	int NodeNum = va_arg(ap, int);

	if (argCnt >= 3)
		Mark = va_arg(ap, int*);

	va_end(ap);

	//cout << "in rrr" << endl;
	if (NodeNum < n / 1000)
	{
		//cout << "in 222" << endl;
		//printf("NodeNum = %d n = %d\n", NodeNum, n);
		sg->n = NodeNum, sg->e = 0;
		int* lab = new int[n], cnt = 0;
		//cout << "in ttt" << endl;
		for (int i = 0; i < sg->n; i++)
		{
			tmploc[cnt] = loc[nodes[i]];
			lab[nodes[i]] = cnt++;
		}
		//cout << "in 666" << endl;
		memcpy(loc, tmploc, sizeof(int) * NodeNum);

		if (argCnt >= 3)
		{
			cnt = 0;
			for (int i = 0; i < sg->n; i++) Mark[cnt++] = nodes[i];
		}

		//cout << "in 999" << endl;

		for (int i = 0; i < NodeNum; i++)
		{
			for (int j = i + 1; j < NodeNum; j++)
			{
				for (int p = cd[nodes[i]]; p < cd[nodes[i] + 1]; p++)
				{
					if (adj[p] == nodes[j])
					{
						sg->e++;
					}
				}

			}
		}
		sg->edges = new edge[sg->e];
		sg->e = 0;

		//cout << "in aaa" << endl;
		for (int i = 0; i < NodeNum; i++)
		{
			for (int j = i + 1; j < NodeNum; j++)
			{
				for (int p = cd[nodes[i]]; p < cd[nodes[i] + 1]; p++)
				{
					if (adj[p] == nodes[j])
					{
						sg->edges[sg->e].s = lab[nodes[i]];
						sg->edges[sg->e].t = lab[nodes[j]];
						sg->e++;
					}
				}

			}
		}
		delete[] lab;
		//cout << "in www" << endl;
		sg->mkGraph();
		//cout << "in eee" << endl;
		return sg;
	}





	sg->n = NodeNum, sg->e = 0;
	int* newFg = new int[n]();

	for (int i = 0; i < NodeNum; i++) newFg[nodes[i]] = 1;

	for (int i = 0; i < e; i++)
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1) sg->e++;

	//printf("sg.e = %d\n", sg.e);
	sg->edges = new edge[sg->e];



	//sort(nodes, nodes+NodeNum);


	/*
	for (int i = 0; i < sg.n; i++)
	{
		int ind = cd[nodes[i]];
		for (int j = i+1; j < sg.n; j++)
		{
			while (ind < cd[nodes[i]] + deg[nodes[i]])
			{
				if (adj[ind] == nodes[j])
				{
					sg.e++;
					ind++;
					break;
				}
				else if (adj[ind] > nodes[j]) break;
				ind++;
			}
		}
	}
	printf("sg.e = %d\n", sg.e);
	sg.edges.resize(sg.e);
	*/

	sg->e = 0;
	int* lab = new int[n], cnt = 0;
	for (int i = 0; i < sg->n; i++)
	{
		tmploc[cnt] = loc[nodes[i]];
		lab[nodes[i]] = cnt++;
	}
	memcpy(loc, tmploc, sizeof(int) * NodeNum);

	if (argCnt >= 3)
	{
		cnt = 0;
		for (int i = 0; i < sg->n; i++) Mark[cnt++] = nodes[i];
	}
	//for (int i = 0; i < sg.n; i++) lab[nodes[i]] = -1;

	/*for (int i = 0; i < sg.n; i++)
	{
		lab[nodes[i]] = cnt++;
		mark[cnt] = node[i];

	}
	*/
	for (int i = 0; i < e; i++)
	{
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
		{
			//lab[edges[i].s] = (lab[edges[i].s] == -1 ? (cnt++) : lab[edges[i].s]);
			//lab[edges[i].t] = (lab[edges[i].t] == -1 ? (cnt++) : lab[edges[i].t]);
			sg->edges[sg->e].s = lab[edges[i].s];
			sg->edges[sg->e].t = lab[edges[i].t];
			sg->e++;
		}
	}


	/*

	for (int i = 0; i < sg.n; i++)
	{
		int ind = cd[nodes[i]];
		for (int j = i + 1; j < sg.n; j++)
		{
			while (ind < cd[nodes[i]] + deg[nodes[i]])
			{
				if (adj[ind] == nodes[j])
				{
					//cout << "nodes[i] = " << nodes[i] << " nodes[j] = " << nodes[j] << endl;
					lab[nodes[i]] = (lab[nodes[i]] == -1 ? (++cnt) : lab[nodes[i]]);
					lab[adj[ind]] = (lab[adj[ind]] == -1 ? (++cnt) : lab[adj[ind]]);
					sg.edges[sg.e].s = lab[nodes[i]];
					sg.edges[sg.e].t = lab[adj[ind]];
					sg.e++;
					ind++;
					break;
				}
				else if (adj[ind] > nodes[j]) break;
				ind++;
			}
		}
	}

	*/
	//printf("sg labeled\n");

	delete[] newFg;
	delete[] lab;
	sg->mkGraph();
	return sg;
}

/*
Graph* Graph::mksubMark(int* nodes, int NodeNum, int *mark)
{
	Graph* sg = new Graph;
	sg->n = NodeNum, sg->e = 0;
	int* newFg = new int[n]();

	for (int i = 0; i < NodeNum; i++) newFg[nodes[i]] = 1;

	for (int i = 0; i < e; i++)
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1) sg->e++;
	sg->edges.resize(sg->e);

	sg->e = 0;
	int* lab = new int[n], cnt = 0;
	//for (int i = 0; i < sg.n; i++) lab[nodes[i]] = -1;

	for (int i = 0; i < sg->n; i++)
	{
		lab[nodes[i]] = cnt;
		mark[cnt++] = nodes[i];
	}

	for (int i = 0; i < e; i++)
	{
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
		{
			//lab[edges[i].s] = (lab[edges[i].s] == -1 ? (cnt++) : lab[edges[i].s]);
			//lab[edges[i].t] = (lab[edges[i].t] == -1 ? (cnt++) : lab[edges[i].t]);
			//mark[lab[edges[i].s]] = edges[i].s;
			//mark[lab[edges[i].t]] = edges[i].t;
			sg->edges[sg->e].s = lab[edges[i].s];
			sg->edges[sg->e].t = lab[edges[i].t];
			sg->e++;
		}
	}

	delete[] newFg;
	delete[] lab;
	sg->mkGraph();
	return sg;
}
*/




int Graph::color(int* color)
{
	iddeg* ig = new iddeg[n];
	for (int i = 0; i < n; i++)
	{
		ig[i].id = i;
		ig[i].degree = deg[i];
	}

	sort(ig, ig + n, IGCmp);

	//color = new int[n];
	memset(color, -1, sizeof(int) * n);
	int* C = new int[(ig[0].degree + 1)]();

	color[ig[0].id] = 0;
	int colorNum = 1;

	for (int i = 1; i < n; i++)
	{
		int tmpDeg = ig[i].degree, tmpid = ig[i].id;
		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < ig[0].degree + 1; j++)
			if (C[j] == 0)
			{
				color[ig[i].id] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}

		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 0;
		}

	}
	//printf("color number = %d\n", colorNum);
	delete[] ig;
	delete[] C;
	return colorNum + 1;
}


long long combination(int n, int m)
{
	if (n < m) return 0;
	if (n == m) return 1;
	long long res = 1;
	for (long long i = n, j = 1; i >= n - m + 1; i--, j++) res *= i, res /= j;
	return res;
}

int outList1(int* list1, int* i, int* j, int s, int s1)
{
	if (*i == s - 1)
		return list1[++ * j];
	if (*j == s1 - 1)
		return list1[++ * i];
	if (list1[*i + 1] < list1[*j + 1])
		return list1[++ * i];
	return list1[++ * j];
}
int merging(int s, int* list1, int s1, int* list2, int s2, int* list3)
{
	int i = 0, j = 0, p = -1, q = s - 1, s3 = 0;
	int x = outList1(list1, &p, &q, s, s1), y = list2[0];
	while (i < s1 && j < s2)
	{
		if (x < y)
		{
			x = outList1(list1, &p, &q, s, s1);
			++i;
			//x = list1[++i];
			continue;
		}
		if (y < x)
		{
			y = list2[++j];
			continue;
		}
		list3[s3++] = x;
		x = outList1(list1, &p, &q, s, s1);
		++i;
		//x = list1[++i];
		y = list2[++j];
	}
	return s3;
}

void nodePass(Graph*& g, int k, int* oldND, int oldN, long long* cnt, bool* live, bool* nType)
{
	//int* ck = new int[k + 1];
	int u, v;


	//if (k == 1000)
		//return;

	//long long* tmpCnt = new long long[g->n]();
	//int* book = new int[g->n];

	//memset(tmpCnt, 0, sizeof(long long) * g->n);

	for (int i = 0; i < oldN; i++)
	{


		u = oldND[i];

		int s, t;



		//ck[0] = u;
		s = t = 0;



		//printf("u = %d g->cd[u] = %d %d\n",u, g->cd[u], g->cd[u + 1]);
		for (int j = g->cd[u]; j < g->cd[u + 1]; j++)
		{

			v = g->adj[j];
			//tmpCnt[v] = 0;
			/*
			uadj[s++] = v;
			for (int p = 0; p < oldN; p++)
			{
				if (oldND[p] == v) s--;
			}
			*/
			if (live[v])
			{
				//uadj[s++] = v;// , book[v] = true;
				if (nType[v] == true) uadj[s++] = v;
				else uadjt[t++] = v;

				//if (v == 7) printf("nType[v]=%d\n",nType[7]);
				/*
				int p;
				for (p = 0; p < oldN; p++)
				{
					if (oldND[p] == v)
					{
						//printf("oldND[i]=%d\n", oldND[p]);
						uadjt[t++] = v;
						break;
					//	if (nType[v] == true)
						//	printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@wrong!\n");
					}
				}
				if (p == oldN)
				{
					uadj[s++] = v;
					//if (nType[v] == false)
					{
					//	printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@wrong false cnt[i] = %d v = %d!\n", cnt[v], v);

					//	return;
					}
				}
				*/




			}
		}

		for (int j = s; j < s + t; j++)
		{
			uadj[j] = uadjt[j - s];
		}


		//sort(uadj, uadj + s + t);


		//for (int j = 0; j < s; j++)
		//{
		//	printf("uadj[%d] = %d\n", j, uadj[j]);
		//}
		//return;
		//printf("s = %d t = %d\n", s, t);

		int left = s + t;
		if (s + t > 2 - 2)
		{
			//ck[0] = u;

			for (int j = 0; j < s; j++)
			{

				v = uadj[j];
				/*
				int flag = 0;
				for (int p = 0; p < oldN; p++)
				{
					if (oldND[p] == v)
					{
						flag = 1;
						break;
					}
				}
				if (flag) continue;
				*/
				//if (book[v] == false)
					//continue;
				//cout << "666 s-j = " << s-j << endl;
				//printf("s = %d j = %d left = %d\n",s,j,left);
				int size = merging(s - j, uadj + j, left, &(g->adj[g->cd[v]]), g->deg[v], Merge);

				if (size < 5 - 2)
				{
					left--;
					continue;
				}


				memcpy(oloc, loc, sizeof(int) * size);
				memcpy(oloc, loc, sizeof(int) * size);
				if (k == 1000)
					continue;

				//cout << "000 size = " << size << " g->n = " << g->n << endl;
				//int* oloc = new int[size];
				//cout << "... merge[0] = " << merge[0] << endl;
				memcpy(oloc, loc, sizeof(int) * size);
				//cout << "size = " << size << endl;
				//cout << "999" << endl;

				Graph* newG = g->mksub(3, Merge, size, mark);
				memcpy(loc, oloc, sizeof(int) * size);
				//delete[] oloc;




				//cout << "888 newG->n = " << newG->n << endl;
				//long long* cntSub = new long long[size]();
				long long nSub = 0;
				if (k - 2 > 1)
					newG->kClique(k - 2, &nSub, cntSub);
				else
				{
					nSub = size;
					for (int p = 0; p < size; p++) cntSub[p] = 1;
				}
				//book[v] = false;
				//for (int p = 0; p < size; p++) book[merge[p]] = false;





				tmpCnt[v] += nSub;

				for (int p = 0; p < size; p++)
				{
					tmpCnt[mark[p]] += cntSub[p];
					cntSub[p] = 0;
				}

				delete newG;
				//delete[] cntSub;


				left--;

				//recursion(kmax, 3, merge, size, g, ck, nck);
			}


			for (int j = g->cd[u]; j < g->cd[u + 1]; j++)
			{
				v = g->adj[j];
				if (live[v]) cnt[v] -= tmpCnt[v], tmpCnt[v] = 0;
			}
		}


		live[u] = false;

		//delete[] uadjt;
		//delete[] uadj;
		//delete[] mark;
		//delete[] merge;
	}
	//delete[] ck;
	//delete[] tmpCnt;
}




bool decpKClique(Graph*& oldG, int k, long long lwb)
{
	Graph* newG;
	bool flag = true;
	double** combin = new double* [k];
	int maxD = oldG->maxDeg;
	for (int i = 0; i < k; i++)
		combin[i] = new double[maxD + 1]();
	for (int i = 0; i <= maxD; i++) combin[0][i] = 1.0;
	for (int i = 0; i < k; i++)	combin[i][i] = 1.0;

	for (int i = 2; i <= maxD; i++)
		for (int j = 1; j < i && j < k; j++)
			combin[j][i] = combin[j][i - 1] + combin[j - 1][i - 1];


	printf("maxD = %d combin[k - 1][maxD] = %lf\n", maxD, combin[k - 1][maxD]);

	double uBound = 0.0, *uBd = new double[oldG->n];
	int  less = 0;
	while (true)
	{


		int* col = new int[oldG->n], newN = 0, oldN = 0;
		int colNum = oldG->color(col);
		int* newND = new int[oldG->n], * oldND = new int[oldG->n];
		//delete[] col;
		int* C = new int[oldG->maxDeg + 1]();

		int t1 = 0, t2 = 0;
		for (int i = 0; i < oldG->n; i++)
		{
			int cv = 0, maxCN = 0, secCN = 0, tmpCN, wCol = -1;
			for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
			{
				if (C[col[oldG->adj[j]]] == 0) cv++;
				C[col[oldG->adj[j]]]++;

				tmpCN = C[col[oldG->adj[j]]];

				if (tmpCN > maxCN)
				{
					if (col[oldG->adj[j]] != wCol)
					{
						secCN = maxCN;
						wCol = col[oldG->adj[j]];
					}
					maxCN = tmpCN;
					continue;
				}

				secCN = secCN > tmpCN ? secCN : tmpCN;

			}


			double uLimit = 0;

			if (cv == 2)
				uLimit = maxCN * combin[k - 2][oldG->deg[i] - maxCN] + combin[k - 1][oldG->deg[i] - maxCN];

			if (cv >= 3)
				uLimit = maxCN * secCN * combin[k - 3][oldG->deg[i] - maxCN - secCN]
				+ (maxCN + secCN) * combin[k - 2][oldG->deg[i] - maxCN - secCN]
				+ combin[k - 1][oldG->deg[i] - maxCN - secCN];
			less = less > (oldG->deg[i] - maxCN - secCN) ? less : (oldG->deg[i] - maxCN - secCN);


			if (uLimit < 0)
				printf("over flow!\n");

			uBound = uBound > uLimit ? uBound : uLimit;

			if (uLimit >= lwb * 1.0) newND[newN++] = i;
			else oldND[oldN++] = i;



			for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
				C[col[oldG->adj[j]]] = 0;
		}


		//printf("t1 = %d t2 = %d\n naive upperBound there newN = %d oldG.n = %d\n", t1, t2, newN, oldG->n);
		if (newN != oldG->n)
			printf("naive upperBound newN = %d oldG->n = %d\n", newN, oldG->n);

		if (newN == oldG->n)
		{
			printf("start dp there newN = %d oldG->n = %d\n", newN, oldG->n);

			newN = 0;
			uBound = 0;

			for (int i = 0; i < oldG->n; i++)
			{
				//printf("i = %d\n", i);
				int colN = 0;
				for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
				{
					if (C[col[oldG->adj[j]]] == 0) colN++;
					C[col[oldG->adj[j]]]++;
				}

				if (colN >= k - 1)
					//if (colN == k - 1) upperBound = 1;
				{
					int* nbrCol = new int[colN];
					colN = 0;

					for (int j = 0; j < oldG->maxDeg + 1; j++)
						if (C[j]) nbrCol[colN++] = C[j];


					double** dpCol = new double* [colN + 1];
					for (int j = 0; j < colN + 1; j++)
						dpCol[j] = new double[k];

					for (int j = 0; j < colN + 1; j++)	dpCol[j][0] = 1;
					for (int j = 1; j < k; j++)			dpCol[j][j] = nbrCol[j - 1] * dpCol[j - 1][j - 1];
					for (int j = 2; j < colN + 1; j++)
						for (int p = 1; p < j && p < k; p++)
							dpCol[j][p] = nbrCol[j - 1] * dpCol[j - 1][p - 1] + dpCol[j - 1][p];


					if (dpCol[colN][k - 1] < 0)
						printf("dpCol over flow\n");

					uBound = uBound > dpCol[colN][k - 1] ? uBound : dpCol[colN][k - 1];
					uBd[i] = dpCol[colN][k - 1];

					if (dpCol[colN][k - 1] >= lwb * 1.0)
						newND[newN++] = i;


					delete[] nbrCol;

					for (int j = 0; j < colN + 1; j++)
						delete[] dpCol[j];
					delete[] dpCol;
				}
				for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
					C[col[oldG->adj[j]]] = 0;

			}

		}

		delete[] C;
		delete[] col;
		if (newN == 0)
		{
			delete[] newND;
			delete oldG;
			flag = false;
			break;
		}
		if (newN < oldG->n)
		{


			newG = oldG->mksub(2, newND, newN);

			delete[] newND;
			delete oldG;
			oldG = newG;
		}
		else
		{
			delete[] newND;
			break;
		}
	}

	//printf("less = %d uBound = %lf\n", less, uBound);

	for (int i = 0; i < k; i++) delete[] combin[i];
	delete[] combin;
	delete[] uBd;

	return flag;
	//Binary(lwb-1, uBound, *oldG, uBd);

}

void Graph::kCliqueCount(int l, long long* tol,
	int* ver, int* lab, int* cdv, int* adjv, int* ns, int** degS, int** subS, long long* cnt)
{
	int u, v, w, end;
	if (l == 2)
	{
		int k = _msize(ver) / sizeof(4) - 1;
		for (int i = 0; i < ns[2]; i++)
		{//list all edges
			u = subS[2][i];
			ver[2] = u;
			//(*n)+=g->d[2][u];
			end = cdv[u] + degS[2][u];
			for (int p = 2; p <= k; p++)
				cnt[ver[p]] += degS[2][u];

			(*tol) += degS[2][u];

			for (int j = cdv[u]; j < end; j++)
			{
				//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
				cnt[adjv[j]]++;
			}

		}
		return;
	}

	for (int i = 0; i < ns[l]; i++)
	{
		u = subS[l][i];
		ver[l] = u;
		//printf("%u %u\n",i,u);
		ns[l - 1] = 0;
		end = cdv[u] + degS[l][u];
		for (int j = cdv[u]; j < end; j++) //relabeling nodes and forming U'.
		{
			v = adjv[j];
			if (lab[v] == l) {
				lab[v] = l - 1;
				subS[l - 1][ns[l - 1]++] = v;
				degS[l - 1][v] = 0;//new degrees
			}
		}
		for (int j = 0; j < ns[l - 1]; j++) //reodering adjacency list and computing new degrees
		{
			v = subS[l - 1][j];
			end = cdv[v] + degS[l][v];
			for (int k = cdv[v]; k < end; k++)
			{
				w = adjv[k];
				if (lab[w] == l - 1)
					degS[l - 1][v]++;
				else
				{
					adjv[k--] = adjv[--end];
					adjv[end] = w;
				}
			}
		}

		kCliqueCount(l - 1, tol, ver, lab, cdv, adjv, ns, degS, subS, cnt);

		for (int j = 0; j < ns[l - 1]; j++) {//restoring labels
			v = subS[l - 1][j];
			lab[v] = l;
		}

	}


}

void Graph::kClique(int k, long long* tol, long long* cnt)
{
	// ord_core
	coreDecomposition();

	// relabel
	int sCore, tCore;
	for (int i = 0; i < e; i++)
	{
		if (coreNum[edges[i].s] > coreNum[edges[i].t])
		{
			edges[i].s = edges[i].s ^ edges[i].t;
			edges[i].t = edges[i].s ^ edges[i].t;
			edges[i].s = edges[i].s ^ edges[i].t;
		}
	}

	// mkspecial
	int* d, * sub, * lab, * cdv, * adjv, * ns, ** degS, ** subS;
	int nsg, maxDv;
	d = new int[n]();
	for (int i = 0; i < e; i++)	d[edges[i].s]++;

	cdv = new int[n + 1];
	nsg = 0, maxDv = 0, cdv[0] = 0;
	sub = new int[n], lab = new int[n];

	for (int i = 1; i < n + 1; i++)
	{
		cdv[i] = cdv[i - 1] + d[i - 1];
		maxDv = (maxDv > d[i - 1]) ? maxDv : d[i - 1];
		sub[nsg++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}

	adjv = new int[e];
	for (int i = 0; i < e; i++)
		adjv[cdv[edges[i].s] + d[edges[i].s]++] = edges[i].t;

	ns = new int[k + 1];
	ns[k] = nsg;

	degS = new int* [k + 1], subS = new int* [k + 1];

	for (int i = 2; i < k; i++)
	{
		degS[i] = new int[n];
		subS[i] = new int[maxDv];
	}
	degS[k] = d;
	subS[k] = sub;

	int* ver = new int[k + 1];
	kCliqueCount(k, tol, ver, lab, cdv, adjv, ns, degS, subS, cnt);


	delete[] lab, delete[] cdv, delete[] adjv, delete[] ns;

	for (int i = 2; i <= k; i++)
		delete[] degS[i], delete[] subS[i];
	delete[] degS, delete[] subS;
}


long long Binary(double left, double right, Graph& g, double* uBd)
{
	//int left = 0, right = n;      //解的范围初始为(0,n],不包含0                             
	while (left + 1 < right)
	{
		long long mid = (right + left) / 2;
		if (1)//check()
		{
			right = mid;          //修正解的范围(left,right]
		}
		else
			left = mid;
	}
	return right;                //最后left + 1 = right
}


bool check(Graph*& oldG, int k, long long lwb, bool* live, long long* cnt, int* liveN, int* rm, int* rmN, int* rmMark)
{
	//Graph* oldG = new Graph(g), 
	Graph* newG;
	int* oldND = new int[oldG->n], lastN = *liveN;
	bool* nType = new bool[oldG->n];
	//for (int i = 0; i < oldG->n; i++) nType[i] = true;

	while (true)
	{

		int* newND, newN = 0, oldN = 0, t = 0;
		//newND = new int[oldG->n];
		lastN = oldN;
		long long sumcnt = 0;
		for (int i = 0; i < oldG->n; i++)
		{

			if (live[i])
			{
				sumcnt += cnt[i];
				t++;
				if (cnt[i] < lwb)
					oldND[oldN++] = i, nType[i] = false, rm[(*rmN)++] = i, rmMark[i] = (*rmN) - 1;
				else nType[i] = true;

			}

		}

		printf("check: n = %d oldN = %d\n", t, oldN);

		if (sumcnt < t * lwb)
		{
			printf("sumcnt check\n");
			return false;
		}



		if (oldN == 0 || oldN == *liveN)
		{
			delete[] oldND;
			//delete oldG;
			if (oldN == 0)
				return true;
			return false;
		}

		*liveN -= oldN;
		nodePass(oldG, k, oldND, oldN, cnt, live, nType);




		/*


		printf("\ngraph -----");
		long long n = 0;
		long long* cnt = new long long[oldG->n]();
		oldG->kClique(k, &n, cnt);

		for (int i = 0; i < oldG->n; i++)
			if (cnt[i] >= lwb) newND[newN++] = i;

*/


/*
for (int i = 0; i < oldG->n; i++)
{
	Graph sg;
	int* iNbr = new int[oldG->deg[i]];
	memcpy(iNbr, oldG->adj + oldG->cd[i], sizeof(int) * oldG->deg[i]);

	oldG->mksub(sg, iNbr, oldG->deg[i]);
	delete[] iNbr;

	long long n = 0;
	sg.kClique(k-1, &n);
	//printf("i = %d n = %lld\n", i, n);


	if (n >= lwb) newND[newN++] = i;
}
*/



/*
delete[] cnt;
if (newN == 0 || newN == oldG->n)
{
	delete[] newND;
	//delete oldG;
	if (newN == 0)
		return false;
	return true;
}




printf("new graph.n = %d", newN);
newG = new Graph;
oldG->mksub(*newG, newND, newN);
delete[] newND;
delete oldG;
oldG = newG;
*/
	}

}



bool decp(Graph*& oldG, int k, long long lwb, long long*& cnt, bool*& live)
{
	//Graph* oldG = new Graph(g), 
	Graph* newG;
	while (true)
	{
		int* newND, newN = 0, *oldND, oldN = 0;
		bool* nType = new bool[oldG->n];
		newND = new int[oldG->n], oldND = new int[oldG->n];

		printf("\ngraph -----");
		long long n = 0;
		//long long* cnt = new long long[oldG->n]();
		oldG->kClique(k, &n, cnt);

		long long maxcnt = 0;
		for (int i = 0; i < oldG->n; i++)
			if (cnt[i] >= lwb) newND[newN++] = i, maxcnt = maxcnt > cnt[i] ? maxcnt : cnt[i], nType[i] = true;
			else oldND[oldN++] = i, nType[i] = false;


		printf("maxcnt = %lld\n", maxcnt);


		/*
		for (int i = 0; i < oldG->n; i++)
		{
			Graph sg;
			int* iNbr = new int[oldG->deg[i]];
			memcpy(iNbr, oldG->adj + oldG->cd[i], sizeof(int) * oldG->deg[i]);

			oldG->mksub(sg, iNbr, oldG->deg[i]);
			delete[] iNbr;

			long long n = 0;
			sg.kClique(k-1, &n);
			//printf("i = %d n = %lld\n", i, n);


			if (n >= lwb) newND[newN++] = i;
		}
		*/
		//delete[] cnt;
		if (newN == 0 || newN == oldG->n)
		{
			delete[] newND;
			//delete oldG;
			if (newN == 0)
				return false;
			return true;
		}
		printf("new graph.n = %d\n", newN);

		if (newN > oldG->n / 2)
		{
			/*
			bool* live = new bool[oldG->n];
			for (int i = 0; i < oldN; i++)
			{
				live[i] = true;
			}
			*/
			nodePass(oldG, k, oldND, oldN, cnt, live, nType);
		}
		else
		{
			newG = oldG->mksub(2, newND, newN);
			//int* mark = new int[newN];

			delete[] newND;
			delete oldG;
			oldG = newG;
			delete[] cnt;
			delete[] live;

			cnt = new long long[oldG->n]();
			live = new bool[oldG->n];
			for (int i = 0; i < oldG->n; i++) live[i] = true;
			long long n = 0;
			oldG->kClique(k, &n, cnt);

		}


		/*
		tmpcnt = new long long[newN];
		for (int i = 0; i < newN; i++)
			tmpcnt[i] = cnt[mark[i]];
		delete[] cnt;
		delete[] mark;
		//return tmpcnt;
		*/

		//decpKClique(oldG, k, lwb);
		break;
	}
}

void kCliqueCore(Graph* dg, int* dgLoc, Graph*& oldG, int k, long long lwb, int N)
{
	decpKClique(oldG, k, lwb);
	bool* live = new bool[oldG->n], * cplive = new bool[oldG->n];
	long long* cnt = new long long[oldG->n]();
	long long* cpcnt = new long long[oldG->n]();
	for (int i = 0; i < oldG->n; i++) live[i] = true;
	decp(oldG, k, lwb, cnt, live);

	long long n = 0;


	//oldG->kClique(k, &n, cnt);

	//for (int i = 0; i < oldG->n; i++) live[i] = true;

	memcpy(cpcnt, cnt, sizeof(long long) * oldG->n);
	memcpy(cplive, live, sizeof(bool) * oldG->n);

	Graph* newG;
	int times = 1, liveN = oldG->n, cpliveN = oldG->n;
	long long left = lwb, right;
	int* tLoc = new int[oldG->n];

	while (true)
	{
		//newG = new Graph(*oldG);
		long long nowb = times * lwb;
		/*
		int wcnt = 0;
		for (int i = 0; i < oldG->n; i++)
		{
			if (combination(oldG->deg[i], k-1) >= nowb)
				wcnt++;
		}
		printf("wcnt = %d\n",wcnt);
		*/


		/*
		if (times == 4)
		{

			if (decpKClique(oldG, k, nowb) == false)
			{
				printf("decpKClique false nowb = %lld\n", nowb);
				//right = nowb;
				//oldG = newG;
				return;
				//break;
			}
		}
		*/
		int* rm = new int[oldG->n], * rmN = new int, * rmMark = new int[oldG->n];
		for (int i = 0; i < oldG->n; i++) rmMark[i] = -1;
		*rmN = 0;
		long long* rmNck = new long long[oldG->n];
		Graph* tmpG;

		if (check(oldG, k, nowb, live, cnt, &liveN, rm, rmN, rmMark))
		{


			if (times == 1)
			{
				ofstream dFile("kcResult.txt", ios::out);
				long long n = 0, *dgCnt = new long long[dg->n](), *rmNck0 = new long long[dg->n];
				dg->kClique(k, &n, dgCnt);


				bool* del = new bool[N], * dgLive = new bool[dg->n], * delTp = new bool[dg->n];
				for (int i = 0; i < N; i++) del[i] = true;
				for (int i = 0; i < oldG->n; i++)
					if (rmMark[i] < 0) del[loc[i]] = false;

				int* rm0 = new int[dg->n], rmN0 = 0, *rmMark0 = new int[dg->n];
				for (int i = 0; i < dg->n; i++)
				{
					dgLive[i] = true;
					if (del[dgLoc[i]] == true)
					{
						delTp[i] = true;
						rm0[rmN0] = i;
						rmNck0[rmN0] = dgCnt[i];
						rmMark0[i] = rmN0;
						rmN0++;
					}
					else delTp[i] = false;
				}

				bheapLLU* heap = mkheapLLU(rmNck0, rmN0);
				keyvalueLLU kv;

				//cout << "111" << endl;

				for (int i = 0; i < rmN0; i++)
				{
					kv = popminLLU(heap);
					//printf("i = %d rmN0 = %d times = 1:\nkv.key = %d rm = %d id = %d kCliqueCoreNumber = %lld\n\n", i, rmN0, kv.key, rm0[kv.key], dgLoc[rm0[kv.key]], kv.value);
					printf("\ri = %d", i);
					//dFile << dgLoc[rm0[kv.key]] << " " << kv.value << endl;

					kCCN = max(kCCN, kv.value);
					cv[cvCnt].key = dgLoc[rm0[kv.key]];
					cv[cvCnt++].value = kCCN;


					//cout << "444" << endl;
					int u = rm0[kv.key];
					//memcpy(tmpc2, tmpc,sizeof(long long) * oldG->n);
					nodePass(dg, 1000, &u, 1, dgCnt, dgLive, delTp);
					//cout << "222" << endl;

					for (int j = dg->cd[u]; j < dg->cd[u + 1]; j++)
					{
						int v = dg->adj[j];
						if (delTp[v] == true && heap->pt[rmMark0[v]] != -1)
						{
							//printf("v = %d rmMark[v] = %d cpcnt[v] - tmpc[v] = %lld\n",v, rmMark[v], cpcnt[v] - tmpc[v] );
							updateLLU(heap, rmMark0[v], dgCnt[v]);
							//updateLLU failed ???
						}
					}
					//cout << "333" << endl;
				}




				/*
				int* ind = new int[N], leftN = dg->n - liveN;

				for (int i = 0; i < dg->n; i++) ind[dgLoc[i]] = i;

				bool* nType = new bool[dg->n];

				*/

				dFile.close();
				return;

			}



			if (times == 2)
			{
				for (int i = 0; i < *rmN; i++) rmNck[i] = cpcnt[rm[i]];
				bool* nType = new bool[oldG->n];
				for (int i = 0; i < oldG->n; i++)
					if (rmMark[i] >= 0) nType[i] = true;
					else nType[i] = false;

				long long* tmpc = new long long[oldG->n];
				memcpy(tmpc, cpcnt, sizeof(long long) * oldG->n);

				bheapLLU* heap = mkheapLLU(rmNck, *rmN);
				keyvalueLLU kv;

				for (int i = 0; i < *rmN; i++)
				{
					kv = popminLLU(heap);
					printf("kv.key = %d rm = %d id = %d kCliqueCoreNumber = %lld\n", kv.key, rm[kv.key], loc[rm[kv.key]], kv.value);

					kCCN = max(kCCN, kv.value);
					cv[cvCnt].key = loc[rm[kv.key]];
					cv[cvCnt++].value = kCCN;



					int u = rm[kv.key];
					//memcpy(tmpc2, tmpc,sizeof(long long) * oldG->n);
					nodePass(oldG, k, &u, 1, tmpc, cplive, nType);


					for (int j = oldG->cd[u]; j < oldG->cd[u + 1]; j++)
					{
						int v = oldG->adj[j];
						if (rmMark[v] >= 0 && heap->pt[rmMark[v]] != -1)
						{
							//printf("v = %d rmMark[v] = %d cpcnt[v] - tmpc[v] = %lld\n",v, rmMark[v], cpcnt[v] - tmpc[v] );
							updateLLU(heap, rmMark[v], tmpc[v]);
							//updateLLU failed ???
						}

					}


				}
			}





			for (int i = 0; i < oldG->n; i++) rmMark[i] = -1;




			printf("check true nowb = %lld\n", nowb);
			left = nowb;
			times *= 2;

			cpliveN = liveN;

			int* nodes = new int[liveN], nodeN = 0, *mark = new int[liveN];
			for (int i = 0; i < oldG->n; i++)
			{
				if (live[i]) nodes[nodeN++] = i;
				else live[i] = true;
			}
			Graph* newG = oldG->mksub(3, nodes, liveN, mark);
			memcpy(tLoc, loc, sizeof(int) * oldG->n);

			long long* tmpcnt = new long long[liveN];
			for (int i = 0; i < liveN; i++)
			{
				tmpcnt[i] = cnt[mark[i]];
			}
			delete[] cnt;
			cnt = tmpcnt;
			delete oldG;
			oldG = newG;
			delete[] nodes;
			delete[] mark;
			memcpy(cplive, live, sizeof(bool) * oldG->n);
			memcpy(cpcnt, cnt, sizeof(long long) * oldG->n);



			//delete newG;
		}
		else
		{

			printf("check false nowb = %lld\n", nowb);
			right = nowb;
			memcpy(live, cplive, sizeof(bool) * oldG->n);
			memcpy(cnt, cpcnt, sizeof(long long) * oldG->n);
			memcpy(loc, tLoc, sizeof(int) * oldG->n);
			liveN = cpliveN;

			//delete oldG;
			//oldG = newG;
			break;
		}
	}


	int maxCoreN = liveN;
	long long maxCore = 0;
	int* pnodes = new int[liveN], pn;
	bool* nType = new bool[oldG->n];
	while (true)
	{
		pn = 0;
		long long minCnt = 0x7fffffffffffffff;// LLONG_MAX;
		for (int i = 0; i < oldG->n; i++)
		{
			nType[i] = true;
			if (live[i] && cnt[i] < minCnt) minCnt = cnt[i];
		}
		//printf("maxCore = %lld minCnt = %lld\n", maxCore, minCnt);
		if (minCnt > maxCore)
		{

			maxCore = minCnt;
			maxCoreN = liveN;
		}
		for (int i = 0; i < oldG->n; i++)
		{
			if (live[i] && cnt[i] == minCnt)
			{
				pnodes[pn++] = i;
				nType[i] = false;
				kCCN = max(kCCN, maxCore);
				cv[cvCnt].key = loc[i];
				cv[cvCnt++].value = kCCN;
				printf("i = %d loc[i] = %d value = %lld\n", i, loc[i], maxCore);
			}
		}

		liveN -= pn;
		//printf("liveN = %d\n", liveN);
		if (liveN == 0)
		{
			break;
		}
		nodePass(oldG, k, pnodes, pn, cnt, live, nType);
		//
	}

	printf("maxCore = %lld maxCoreN = %d\n", maxCore, maxCoreN);

	/*
	while (left + 1 < right)
	{
		newG = new Graph(*oldG);
		long long mid = (right + left) / 2;
		//decpKClique(oldG, k, mid);
		if (check(oldG, k, mid, live, cnt, &liveN))//check()
		{
			printf("check true mid = %lld\n", mid);
			left = mid;          //修正解的范围[left,right)
			//memcpy(cplive, live, sizeof(bool) * oldG->n);
			//memcpy(cpcnt, cnt, sizeof(long long) * oldG->n);
			cpliveN = liveN;

			int* nodes = new int[liveN], nodeN = 0, * mark = new int[liveN];
			for (int i = 0; i < oldG->n; i++)
			{
				if (live[i]) nodes[nodeN++] = i;
				else live[i] = true;
			}
			Graph* newG = new Graph;
			oldG->mksubMark(*newG, nodes, liveN, mark);
			long long* tmpcnt = new long long[liveN];
			for (int i = 0; i < liveN; i++)
			{
				tmpcnt[i] = cnt[mark[i]];
			}
			delete[] cnt;
			cnt = tmpcnt;
			delete oldG;
			oldG = newG;
			delete[] nodes;
			delete[] mark;
			memcpy(cplive, live, sizeof(bool) * oldG->n);
			memcpy(cpcnt, cnt, sizeof(long long) * oldG->n);


			//delete newG;
		}
		else
		{
			printf("check false mid = %lld\n", mid);
			right = mid;
			memcpy(live, cplive, sizeof(bool) * oldG->n);
			memcpy(cnt, cpcnt, sizeof(long long) * oldG->n);
			liveN = cpliveN;
			//delete oldG;
			//oldG = newG;
		}
	}

	*/

	/*
	int num = 0;
	long long aveCoreNum = 0;
	for (int i = 0; i < oldG->n; i++)
	{
		if (live[i])
		{
			num++;
			aveCoreNum += cnt[i];
		}

	}
	aveCoreNum /= k;
	printf("End: n = %d aveCoreNum = %lld\n", num, aveCoreNum/num);
	*/

	//return left;

}


int calDP(Graph& g, int i)
{

}

int main(int argc, char** argv)
{
	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	Graph g;
	int k = atoi(argv[1]);
	cout << "Reading edgelist from file " << argv[2] << endl;
	g.readedgelist(argv[2]);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;
	int* color = new int[g.n];

	int colorNum = g.color(color);



	double** dp = new double* [g.n];


	//double*** ck = new double** [g.n];

	double** ck = new double* [colorNum + 1];
	for (int j = 0; j < colorNum + 1; j++)
	{
		ck[j] = new double[k + 1];
	}
	printf("colorNum = %d\n", colorNum);

	int** CC = new int* [g.n];

	//int** nbrCol = new int* [g.n];

	for (int i = 0; i < g.n; i++)
	{
		dp[i] = new double[k + 1];
		int colorNum_i = 0;
		//int* C = new int[colorNum]();
		CC[i] = new int[colorNum]();
		for (int j = g.cd[i]; j < g.cd[i + 1]; j++)
		{
			int nbr = g.adj[j];
			CC[i][color[nbr]]++;
			colorNum_i = max(colorNum_i, color[nbr]);
		}
		colorNum_i++;
		double* NotColor0 = new double[k + 1]();
		double* MustColor0 = new double[k + 1]();
		//int* tol = new int[k + 1];

		NotColor0[0] = 1;
		for (int c = 1; c < colorNum_i; c++)
		{
			for (int j = k; j > 0; j--)
				NotColor0[j] = NotColor0[j - 1] * CC[i][c] + NotColor0[j];
		}

		for (int j = 1; j <= k; j++)
		{
			MustColor0[j] = NotColor0[j - 1] * CC[i][0];
			dp[i][j] = MustColor0[j] + NotColor0[j];
		}
		delete[] NotColor0;
		delete[] MustColor0;



		//------------------------------------------------------------------

		int tolCol = colorNum;
		//double** ck = new double* [tolCol + 1];
		for (int j = 0; j < tolCol + 1; j++)
		{
			//ck[j] = new double[k + 1];
			ck[j][0] = 1;
		}
		for (int j = 0; j < k + 1; j++)
			ck[0][j] = 0;

		ck[0][0] = 1;


		//nbrCol[i] = new int[tolCol]();

		for (int j = g.cd[i]; j < g.cd[i + 1]; j++)
		{
			int v = g.adj[j];
			//nbrCol[i][color[v]]++;
		}


		int cNum = tolCol;


		for (int j = 1; j < cNum + 1; j++)
			for (int p = 1; p < k + 1; p++)
				ck[j][p] = ck[j - 1][p - 1] * CC[i][j - 1] + ck[j - 1][p];


		//delete[] C;

		//if (i == 2660)
		//{
		//	printf("%d %d \n", cNum, colorNum_i);
		//	for (int q = 0; q < colorNum_i; q++)
		//	{
		//		printf("%d ", C[q]);

		//	}
		//	printf("\n");
		//	for (int q = 0; q < colorNum_i; q++)
		//	{
		//		printf("%d ", nbrCol[i][q]);

		//	}
		//	printf("\n");
		if (dp[i][k] < 0)
			printf("%d %lf %lf\n", i, dp[i][k], ck[cNum][k]);

		//if (dp[i][k] != ck[i][cNum][k])
			//printf("Not %d colorNum = %d %lf %lf\n", i, colorNum_i, dp[i][k], ck[i][cNum][k]);

		//}
	}

	printf("times 111\n");




	int delNodes = 0, times = 0;
	bool* mark = new bool[g.n]();
	int leftN = g.n, leftM = g.e;

	double* NotColor0 = new double[k + 1]();
	double* MustColor0 = new double[k + 1]();
	//int* C = new int[colorNum]();

	double tolMax = 0;
	int maxN = 0, maxM = 0;

	while (leftN > 0)
	{
		times++;
		leftN -= delNodes;
		if (leftN < 0) printf("%d %d\n", leftN, delNodes);
		if (leftN <= 0) break;

		delNodes = 0;
		double Min = 1e300;
		for (int i = 0; i < g.n; i++)
		{
			if (mark[i] == true) continue;
			Min = min(Min, dp[i][k]);
		}
		if (Min > tolMax)
		{
			tolMax = Min;
			maxN = leftN;
			maxM = leftM;
			//printf("tolMax = %lf maxN = %d\n", tolMax, maxN);
		}
		//tolMax = max(tolMax, Min);

		if (times % 500 == 0)
			printf("times = %d left nodes = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", times, leftN, tolMax, maxN, maxM, 1.0 * maxM / maxN);


		for (int i = 0; i < g.n; i++)
		{
			if (mark[i] == true) continue;
			if (fabs(dp[i][k] - Min) < 1e-9)
			{
				mark[i] = true;
				delNodes++;
				leftM -= g.deg[i];

				//printf("delNodes = %d\n", delNodes);

				for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
				{
					int nbr = g.adj[j];

					for (int h = g.cd[nbr]; h < g.cd[nbr] + g.deg[nbr]; h++)
					{
						int hnbr = g.adj[h];
						if (hnbr == i)
						{
							swap(g.adj[h], g.adj[g.cd[nbr] + g.deg[nbr] - 1]);
							g.deg[nbr]--;
							break;
						}
					}





					//for (int ii = 0; ii < colorNum; ii++)
					//{
					//	C[ii] = 0;
					//}

					//int colorNum_i = 0;
					//for (int h = g.cd[nbr]; h < g.cd[nbr] + g.deg[nbr]; h++)
					//{
					//	int hnbr = g.adj[h];
					//	CC[nbr][color[hnbr]]++;
					//	colorNum_i = max(colorNum_i, color[hnbr]);
					//}
					//colorNum_i++;

					CC[nbr][color[i]]--;


					//printf("Min = %lf dp = %lf %d \n", Min, dp[i][k], colorNum_i);




					double lastdp = dp[nbr][k];
					NotColor0[0] = 1;
					for (int h = 1; h <= k; h++)
					{
						MustColor0[h] = NotColor0[h - 1] * (CC[nbr][color[i]] + 1);

						NotColor0[h] = dp[nbr][h] - MustColor0[h];
						//if (nbr == 137729)
							//printf("%d %lf %lf\n", color[i], MustColor0[h], NotColor0[h] );
					}

					for (int h = 1; h <= k; h++)
					{
						MustColor0[h] = NotColor0[h - 1] * CC[nbr][color[i]];
						dp[nbr][h] = NotColor0[h] + MustColor0[h];
					}


					//----------------------------------


					/*
					for (int j = 0; j < colorNum + 1; j++)
					{
						ck[j][0] = 1;
					}
					for (int j = 0; j < k + 1; j++)
						ck[0][j] = 0;

					ck[0][0] = 1;

					int cNum = colorNum;


					for (int h = 1; h < cNum + 1; h++)
						for (int p = 1; p < k + 1; p++)
							ck[h][p] = ck[h - 1][p - 1] * C[h - 1] + ck[h - 1][p];
					*/

					//if (nbr == 1)
					//{
					//	for (int q = 0; q < colorNum_i; q++)
					//	{			
					//		printf("%d ", C[q]);

					//	}
					//	printf("\n");
					//	printf("color[i] = %d %lf AAA %d colorNum = %d %lf %lf\n", color[i], lastdp, nbr, colorNum_i, dp[nbr][k], ck[nbr][cNum][k]);
					//}

					//printf("color[i] = %d %lf AAA %d colorNum = %d %lf %lf\n", color[i], lastdp, nbr, colorNum_i, dp[nbr][k], ck[nbr][cNum][k]);

					//if ( ck[nbr][cNum][k] != dp[nbr][k])


					//----------------------------------
						//if(fabs(ck[cNum][k] - dp[nbr][k])/ck[cNum][k] > 1e-9)
							//printf("color[i] = %d %lf AAA %d colorNum = %d %lf %lf\n", color[i], lastdp, nbr, colorNum_i, dp[nbr][k], ck[cNum][k]);



				}
			}
		}

	}
	printf("End: times = %d left nodes = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", times, leftN, tolMax, maxN, maxM, 1.0 * maxM / maxN);


	delete[] NotColor0;
	delete[] MustColor0;
	//delete[] C;






	printf("The End\n");

	return 0;

}
