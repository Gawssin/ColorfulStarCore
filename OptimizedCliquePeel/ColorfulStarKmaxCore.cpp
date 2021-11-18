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
	g.readedgelist(argv2);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;
	auto t1 = getTime();

	int* color = new int[g.n];
	int colorNum = g.color(color);
	printf("colorNum = %d\n", colorNum);

	double** dp = new double* [g.n];
	int** CC = new int* [g.n];
	initColStarDegree(g, dp, h, colorNum, color, CC, 0);

	auto t2 = getTime();

	double maxCore = 0;
	int maxCoreNum = 0;
	double* ColofulStarCoreNum = new double[g.n];
	ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, 0, &maxCore, &maxCoreNum);

	int* maxCoreNodes = new int[maxCoreNum];
	maxCoreNum = 0;
	for (int i = 0; i < g.n; i++)
	{
		if (ColofulStarCoreNum[i] >= maxCore)
			maxCoreNodes[maxCoreNum++] = i;
	}
	Graph& maxCoreSub = g.mksub(2, maxCoreNodes, maxCoreNum);

	auto t4 = getTime();
	printf("startPeeling = %lfs\n", ((double)timeGap(t1, t4)) / 1e6);

	hCliquePeeling(maxCoreSub, h);

	printf("\nColorful %d-star Kmax Core\nNodes:\t\t%d\nEdges:\t\t%d\nKmax:\t\t%lf\nClique-Density:\t%lf\n\n", h, maxCoreSub.n, maxCoreSub.e, maxCore, 1.0 * (maxCoreSub.clique->cliqueNumber) / maxCoreSub.n);

	auto t5 = getTime();

	printf("- Overall time = %lfs\n", ((double)timeGap(t1, t5)) / 1e6);

	return 0;
}
