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
#include <unordered_map>
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
	int algoFlag = 0;
	if (argc >= 4)
	{
		if (strcmp(argv[3], "MaxCore") == 0) algoFlag = 0;
		else if (strcmp(argv[3], "MaxCorePeel") == 0) algoFlag = 1;
	}

	auto t0 = getTime();

	Graph g;
	int h = atoi(argv1);
	cout << "Reading edgelist from file " << argv2 << endl;
	g.readedgelist(argv2);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;
	auto t1 = getTime();

	int* color = new int[g.n];
	int colorNum = g.color(color);
	printf("colorNum = %d\n", colorNum);

	__int128** dp = new __int128* [g.n];
	int** CC = new int* [g.n];
	initColStarDegree(g, dp, h, colorNum, color, CC, 0);

	auto t2 = getTime();

	__int128 maxCore = 0;
	int maxCoreNum = 0;
	__int128* ColofulStarCoreNum = new __int128[g.n];
	ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, 0, &maxCore, &maxCoreNum);

	auto t3 = getTime();
	printf("Colorful Star Core Decomposition Finished. Time = %lfs\n", ((double)timeGap(t1, t3)) / 1e6);

	int* maxCoreNodes = new int[maxCoreNum];
	maxCoreNum = 0;
	for (int i = 0; i < g.n; i++)
	{
		if (ColofulStarCoreNum[i] >= maxCore)
			maxCoreNodes[maxCoreNum++] = i;
	}
	Graph& maxCoreSub = g.mksub(2, maxCoreNodes, maxCoreNum);

	if (algoFlag == 0)
	{
		maxCoreSub.clique = new Clique(maxCoreSub.n, maxCoreSub.e, h);
		long long tol = 0;
		long long* cnt = new long long[maxCoreSub.n]();

		for (int i = 0; i < maxCoreNum; i++) maxCoreNodes[i] = i;

		maxCoreSub.kCliqueNew(h, &tol, cnt, maxCoreNodes, maxCoreNum);

		printf("\nColorful %d-star Kmax Core\nNodes:\t\t%d\nEdges:\t\t%d\nKmax:\t\t%s\nClique-Density:\t%lf\n\n", h, maxCoreSub.n, maxCoreSub.e, _int128_to_str(maxCore), 1.0 * tol / maxCoreSub.n);
		delete[] cnt;
	}
	else
	{
		auto t4 = getTime();
		printf("startPeeling = %lfs\n", ((double)timeGap(t1, t4)) / 1e6);

		hCliquePeeling(maxCoreSub, h);
	}

	auto t5 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t1, t5)) / 1e6);
	return 0;
}
