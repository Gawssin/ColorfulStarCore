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

#define debug 0

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
	double maxCore = 0;
	int maxCoreNum = 0;
	double* ColofulStarCoreNum = new double[g.n];
	ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, &maxCore, &maxCoreNum);


	if (debug)
	{
		for (int i = 0; i < g.n; i++)
		{
			printf("id = %d starDegree = %lf\n", i, ColofulStarCoreNum[i]);
		}
	}


	int* maxCoreNodes = new int[maxCoreNum];
	maxCoreNum = 0;
	for (int i = 0; i < g.n; i++)
	{
		if (ColofulStarCoreNum[i] >= maxCore)
			maxCoreNodes[maxCoreNum++] = i;
	}
	Graph maxCoreSub = g.mksub(2, maxCoreNodes, maxCoreNum);
	maxCoreSub.clique = new Clique(maxCoreSub.n, maxCoreSub.e, h);
	long long tol = 0;
	long long* cnt = new long long[maxCoreSub.n]();
	//g.kClique(h, &tol, cnt);

	for (int i = 0; i < maxCoreNum; i++) maxCoreNodes[i] = i;


	maxCoreSub.kCliqueNew(h, &tol, cnt, maxCoreNodes, maxCoreNum);
	printf("The %d-clique density of the colorful %d-star maxK core: %lf\n", h, h, 1.0 * tol / maxCoreNum);


	auto t3 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t2, t3)) / 1e6);

	printf("The End\n");

	return 0;

}
