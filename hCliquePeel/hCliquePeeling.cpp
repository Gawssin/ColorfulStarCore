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


	//g.coreDecomposition();
	//int largeCliqueSize = g.outLargeClique();
	//printf("largeCliqueSize: %d\n", largeCliqueSize);

	hCliquePeeling(g, h);


	auto tClique = getTime();

	printf("- Overall time = %lfs\n", ((double)timeGap(t1, tClique)) / 1e6);

	return 0;

}
