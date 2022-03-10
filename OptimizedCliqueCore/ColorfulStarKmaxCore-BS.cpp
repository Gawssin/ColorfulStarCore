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

	__int128 lowerBound = dp[0][h - 1], upperBound = 0;

	for (int i = 0; i < g.n; i++)
	{
		lowerBound = min(lowerBound, dp[i][h - 1]);
		upperBound = max(upperBound, dp[i][h - 1]);
		memset(CC[i], 0, colorNum * sizeof(int));
	}
	printf("upperBound: %s\n", _int128_to_str(upperBound));

	printf("lowerBound: %s\n", _int128_to_str(lowerBound));

	__int128 maxColorfulCoreNumber = BinaryMaxCore(g, h, lowerBound, upperBound + 1, CC, dp, color);

	//g.deg[u] > 0 if u in the colorful h-star Kmax core

	printf("maxColorfulCoreNumber: %s\n", _int128_to_str(maxColorfulCoreNumber));

	auto t3 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t1, t3)) / 1e6);
	return 0;
}
