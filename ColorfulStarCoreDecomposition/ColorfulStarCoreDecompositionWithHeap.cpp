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


int main(int argc, char** argv)
{
	char* argv1, * argv2, * argv3;
	argv1 = argv[1], argv2 = argv[2];
	int basicFlag = 0;
	if (argc >= 4 && strcmp(argv[3], "basic") == 0)
	{
		basicFlag = 1;
	}
	cout << "Version: " << ( basicFlag ? "Basic (HStarDP)":"Advanced (HStarCD)" ) << endl << endl;
	
	auto t0 = getTime();

	Graph g;
	int h = atoi(argv1);
	cout << "Reading edgelist from file: " << argv2 << endl;
	g.readedgelist(argv2);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;
	auto t1 = getTime();

	int* color = new int[g.n];
	int colorNum = g.color(color);
	printf("Total number of colors: %d\n", colorNum);

	double** dp = new double* [g.n];
	int** CC = new int* [g.n];
	double* ColofulStarCoreNum = new double[g.n];

	initColStarDegree(g, dp, h, colorNum, color, CC, basicFlag);
	ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, basicFlag);


	auto t2 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t1, t2)) / 1e6);

	return 0;
}
