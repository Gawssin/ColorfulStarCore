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
	argv1 = argv[1], argv2 = argv[2];


	//--------- readCMD begin
	//readCMD rCMD(argc, argv);
	//argv1 = rCMD.read();
	//argv2 = rCMD.read();
	//--------- readCMD end


	auto t0 = getTime();

	Graph g;
	int h = atoi(argv1);
	//h = 3;
	cout << "h = " << h << endl;
	cout << "Reading edgelist from file " << argv2 << endl;

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
	//int delNum = 0, delEdges = 0;
	for (int i = 0; i < g.n; i++)
	{
		if (g.coreNum[i] < largeCliqueSize - 1) 
		{
			//delEdges += g.deg[i];
			//deleteNodes(&g, &i, 1);
			//delNum++;
			delArr[delNum++] = i;
			delDeg += g.deg[i];
		}
	}
	interEdges = g.deleteNodes(delArr, delNum);
	delete[] delArr;
	int delWEdges = delDeg - interEdges / 2;
	printf("Get w core, delNum: %d, delEdges = %d\n", delNum, delWEdges);


	
	int* color = new int[g.n];
	int colorNum = g.color(color);
	printf("colorNum = %d\n", colorNum);


	double** dp = new double* [g.n];
	int** CC = new int* [g.n];
	initColStarDegree(g, dp, h, colorNum, color, CC);

	delNum = 0, delEdges = 0;;
	long long LB = combination(largeCliqueSize - 1, h - 1);
	ColorfulStarCore(g, dp, h, color, CC, (double)LB, delNum, delEdges);
	

	printf("Get ColorfulStar core, delNum: %d, delEdges = %d\n", delNum, delEdges);

	printf("Get ColorfulStar core, N: %d, M: %d\n", g.n - delNum, g.e - delWEdges - delEdges);

	//int nowN = 0, nowE = 0;
	//for (int i = 0; i < g.n; i++)
	//{
	//	if (g.deg[i] != 0 )
	//	{
	//		nowN++;
	//		nowE += g.deg[i];
	//	}
	//}
	//printf("Get ColorfulStar core, N: %d, M: %d\n", nowN, nowE);



	hCliquePeeling(g, h);


	auto tClique = getTime();

	printf("- Overall time = %lfs\n", ((double)timeGap(t1, tClique)) / 1e6);

	return 0;

}
