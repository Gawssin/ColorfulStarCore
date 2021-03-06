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
#include "../header/Graph.hpp"
#include "../header/tool.hpp"
#include "../header/ColorfulStarCore.hpp"
#include "../header/hCliquePeel.hpp"

using namespace std;
int main(int argc, char** argv)
{
	char* argv1, * argv2;
	argv1 = argv[1], argv2 = argv[2];

	auto t0 = getTime();

	Graph g;
	int h = atoi(argv1);

	cout << "Reading edgelist from file " << argv2 << endl;
	cout << "h = " << h << endl;
	g.readedgelist(argv2);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;

	printf("N: %d\tM:%d\n", g.n, g.e);

	auto t1 = getTime();

	g.coreDecomposition();
	int largeCliqueSize = g.outLargeClique();
	printf("largeCliqueSize: %d\n", largeCliqueSize);

	int maxCore = -1;
	for (int i = 0; i < g.n; i++)
	{
		maxCore = max(maxCore, g.coreNum[i]);
	}
	printf("maxCore: %d\n", maxCore);

	int* delArr = new int[g.n], delNum = 0, delEdges = 0, delDeg = 0;
	int interEdges;	//the number of edges connecting u, v which both are in delArr.

	for (int i = 0; i < g.n; i++)
	{
		if (g.coreNum[i] < largeCliqueSize - 1)
		{
			delArr[delNum++] = i;
			delDeg += g.deg[i];
		}
	}
	interEdges = g.deleteNodes(delArr, delNum);
	delete[] delArr;
	int delWEdges = delDeg - interEdges / 2;
	int* color = new int[g.n];
	int colorNum = g.color(color);
	printf("colorNum = %d\n", colorNum);

	__int128** dp = new __int128* [g.n];
	int** CC = new int* [g.n];
	initColStarDegree(g, dp, h, colorNum, color, CC, 0);

	delNum = 0, delEdges = 0;;
	long long LB = combination(largeCliqueSize - 1, h - 1);
	ColorfulStarCore(g, dp, h, color, CC, (__int128)LB, delNum, delEdges);

	printf("\nGet ColorfulStar core\nDeleted Nodes:\t%d\nDeleted Edges:\t%d\nLeft Nodes:\t%d\nLeft Edges:\t%d\n\n",
		delNum, delEdges, g.n - delNum, g.e - delWEdges - delEdges);

	auto reductionTime = getTime();
	printf("- Reduction time = %lfs\n", ((double)timeGap(t1, reductionTime)) / 1e6);

	hCliquePeeling(g, h);

	auto tClique = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t1, tClique)) / 1e6);
	return 0;
}
