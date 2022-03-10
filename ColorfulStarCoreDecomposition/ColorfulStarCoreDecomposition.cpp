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

using namespace std;


int main(int argc, char** argv)
{
	char* argv1, * argv2, * argv3;
	argv1 = argv[1], argv2 = argv[2];
	int basicFlag = 0, colorAlgo = atoi(argv[3]);
	if (argc >= 5 && strcmp(argv[4], "basic") == 0)
	{
		basicFlag = 1;
	}
	cout << "Version: " << (basicFlag ? "Basic (HStarDP)" : "Advanced (HStarCD)") << endl << endl;

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
	int colorNum = g.color(color, colorAlgo);
	printf("Total number of colors: %d\n", colorNum);

	__int128** dp = new __int128* [g.n];
	int** CC = new int* [g.n];
	__int128* ColofulStarCoreNum = new __int128[g.n];

	initColStarDegree(g, dp, h, colorNum, color, CC, basicFlag);
	ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, basicFlag);

	auto t2 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t1, t2)) / 1e6);


	// study nodes' distribution

	// if (argc >= 5)
	// {
	//     unordered_map<__int128, int> umap;
	//     for (int i = 0; i < g.n; i++)
	//     {
	//         umap[ColofulStarCoreNum[i]]++;
	//     }

	//     ofstream outfile(argv[5], ios::out);

	//     for (auto x : umap)
	//         outfile << _int128_to_str(x.first) << " " << x.second << endl;

	//     outfile.close();
	// }

	return 0;
}
