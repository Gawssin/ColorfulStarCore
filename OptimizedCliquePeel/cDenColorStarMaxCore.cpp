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
#include "../header/heapLLU.h"
#include "../header/GraphNew.hpp"
#include "../header/tool.hpp"
#include "../header/ColorfulStarCore.hpp"
using namespace std;

#define debug 0

void deleteNodes(Graph* g, int* delArray, int size)
{

	if (size == 1)
	{
		int u = delArray[0];
		for (int j = g->cd[u]; j < g->cd[u] + g->deg[u]; j++)
		{
			int v = g->adj[j];
			for (int k = g->cd[v]; k < g->cd[v] + g->deg[v]; k++)
			{
				int w = g->adj[k];
				if (w == u)
				{
					g->adj[k] = g->adj[g->cd[v] + g->deg[v] - 1];
					g->adj[g->cd[v] + g->deg[v] - 1] = w;
					g->deg[v]--;
					break;
				}
			}
		}
		g->deg[u] = 0;
		return;
	}
}

int main(int argc, char** argv)
{

	char* argv1, * argv2;
	//argv1 = argv[1], argv2 = argv[2];


	//--------- readCMD begin
	readCMD rCMD(argc, argv);
	argv1 = rCMD.read();
	argv2 = rCMD.read();
	//--------- readCMD end


	auto t0 = getTime();

	Graph g;
	int h = atoi(argv1);
	cout << "Reading edgelist from file " << argv2 << endl;
	g.readedgelist(argv2);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;
	auto t1 = getTime();

	//long long tol;
	//long long * cnt = new long long[g.n];
	//g.kClique(k, &tol, cnt);

	//printf("%d-clique: %lld\n", k, tol);

	//return 0;

	int* color = new int[g.n];
	int colorNum = g.color(color);
	printf("colorNum = %d\n", colorNum);

	if (debug)
	{
		//int setColor[] = { 0,3,0,1,2,0,4,3,1,0 };
		//color = setColor;
		for (int i = 0; i < g.n; i++)
		{
			printf("id = %d color = %d\n", i, color[i]);
		}
	}

	double** dp = new double* [g.n];
	int** CC = new int* [g.n];
	initColStarDegree(g, dp, h, colorNum, color, CC);

	auto t2 = getTime();


	if (debug)
	{
		for (int i = 0; i < g.n; i++)
		{
			printf("id = %d starDegree = %lf\n", i, dp[i][h - 1]);
		}
	}

	double* ColofulStarCoreNum = new double[g.n];
	ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum);

	if (debug)
	{
		for (int i = 0; i < g.n; i++)
		{
			printf("id = %d starDegree = %lf\n", i, ColofulStarCoreNum[i]);
		}
	}

	auto t3 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t2, t3)) / 1e6);

	printf("The End\n");






	//----------------------------------
	int* GNodes = new int[g.n];
	for (int i = 0; i < g.n; i++) GNodes[i] = i;
	g.clique = new Clique(g.n, g.e, h);


	long long tol = 0;
	long long* cnt = new long long[g.n]();
	//g.kClique(h, &tol, cnt);

	g.kCliqueNew(h, &tol, cnt, GNodes, g.n);

	printf("Toltal cliques: %lld\th-clique density of colorful h-star maxK core: %lf\n", tol, tol / 162);



	bheapLLU<long long>* cHeap = mkheapLLU<long long>(g.n, cnt);

	keyvalueLLU<long long> kv;

	long long curCliqueCore = 0, leftClique = tol;
	int leftN = g.n, leftM = g.e;
	int maxCliDenN = 0;	//the number of nodes in the subgraph achiving the largest hclique density;
	int maxCliDenM = 0;	//the number of edges in the subgraph achiving the largest hclique density;

	int maxCliCoreDenN = 0;	//the number of nodes in the maximal hclique core;
	int maxCliCoreDenM = 0;	//the number of edges in the maximal hclique core;
	double cliqueDensity = 0.0, curCliqueDensity, cliqueCoreDen = 0.0;

	int* nbrArr = new int[g.maxDeg], nbrNum;
	long long* nbrCnt = new long long[g.n](), nbrTol;
	long long maxCliDeg = 0;

	int coreTocore = 0, cntctc = 0;

	//peeling ordering infects maxCliqueDensity.

	while (leftN > 0)
	{
		if (leftN == 221)
		{
			printf("leftN = %d, leftClique = %d\n", leftN, leftClique);
		}

		kv = popminLLU<long long>(cHeap);
		long long cliqueDeg = kv.value;
		curCliqueDensity = 1.0 * leftClique / leftN;

		if (maxCliDeg < cliqueDeg)
		{
			cliqueCoreDen = curCliqueDensity;
			maxCliCoreDenN = leftN;
			maxCliCoreDenM = leftM;

			maxCliDeg = cliqueDeg;
			cntctc += coreTocore;
			//printf("maxCliDeg = %lld -> %lld, coreTocore = %d, \t leftN = %d, \t tolcn = %d\n", maxCliDeg - 1, maxCliDeg, coreTocore, leftN, cntctc);

			coreTocore = 0;
		}
		coreTocore++;

		maxCliDeg = max(maxCliDeg, cliqueDeg);

		//printf("id = %d, cliqueDeg = %lld, leftClique = %lld, leftN = %d\n", kv.key, kv.value, leftClique, leftN);



		//if (g.n - leftN < 100)
			//printf("-------------------------maxCliDeg = %lld, delId = %d, leftN = %d, leftClique: %lld\n", maxCliDeg, kv.key, leftN, leftClique);


		if (cliqueDensity < curCliqueDensity)
		{
			cliqueDensity = curCliqueDensity;
			maxCliDenN = leftN;
			maxCliDenM = leftM;
			//cliqueDensity = 1.0 * leftClique / leftN;
		}

		int delId = kv.key;
		nbrNum = 0;
		for (int i = g.cd[delId]; i < g.cd[delId] + g.deg[delId]; i++)
		{
			int adj = g.adj[i];
			nbrArr[nbrNum++] = adj;
		}
		nbrTol = 0;
		g.kCliqueNew(h - 1, &nbrTol, nbrCnt, nbrArr, nbrNum);


		leftClique -= nbrTol;



		for (int i = 0; i < nbrNum; i++)
		{
			int nbrId = nbrArr[i];
			cnt[nbrId] -= nbrCnt[nbrId];
			updateLLU(cHeap, nbrId, cnt[nbrId]);
			//printf("nbrId = %d, nbrCnt = %lld\n",nbrId, nbrCnt[nbrId]);
		}
		//printf("------------\n");


		for (int i = 0; i < nbrNum; i++)
			nbrCnt[nbrArr[i]] = 0;


		g.deleteNodes(&delId, 1);
		leftN--;
	}



	cntctc += coreTocore;
	printf("after coreTocore = %d, tolcn = %d\n", coreTocore, cntctc);

	printf("maxCliDenN = %d, cliqueDensity = %lf\n", maxCliDenN, cliqueDensity);
	printf("maxCliDeg = %lld, maxCliCoreDenN = %d, cliqueCoreDensity = %lf\n", maxCliDeg, maxCliCoreDenN, cliqueCoreDen);


	auto tClique = getTime();

	printf("- Overall time = %lfs\n", ((double)timeGap(t1, tClique)) / 1e6);

	return 0;


	return 0;

}
