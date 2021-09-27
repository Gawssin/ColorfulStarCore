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




	//char *s = strrchr(argv2, '\\');

	char datasetsName[100];

	for (int i = strlen(argv2) - 1; i >= 0; i--)
	{
		if (argv2[i] == '\\')
		{
			for (int j = 0; j < 100; j++)
			{
				if (argv2[j + i + 1] != '.')
					datasetsName[j] = argv2[j + i + 1];
				else 
				{
					datasetsName[j] = '\0';
					break;
				}
			}
			break;
		}
	}

	//printf("%s\n", datasetsName);

	char path[200] = "F:\\datasets\\datasets-n-m\\datasets-scale\\";
	char newPathName[200];

	srand(unsigned(time(0)));

	int onEdges = 0;
	if (onEdges)
	{
		int* newID = new int[g.n];
		for (int i = 0; i < 4; i++)
		{
			memset(newID, -1, sizeof(int) * g.n);
			random_shuffle(g.edges, g.edges + g.e);
			int newE = 0.2 * (i + 1) * g.e;
			int nowID = 0;
			for (int j = 0; j < newE; j++)
			{
				if (newID[g.edges[j].s] == -1) newID[g.edges[j].s] = nowID++;
				if (newID[g.edges[j].t] == -1) newID[g.edges[j].t] = nowID++;
			}

			FILE* fp;
			sprintf(newPathName, "%sEdges\\%s-0.%d-Edge.txt", path, datasetsName, (i + 1) * 2);

			//printf("%s\n", newPathName);

			fp = fopen(newPathName, "w");

			fprintf(fp, "%d %d\n", nowID, newE);

			for (int j = 0; j < newE; j++)
			{
				fprintf(fp, "%d %d\n", newID[g.edges[j].s], newID[g.edges[j].t]);
			}
			fclose(fp);
			printf("%lf onEdges End\n", 0.2 * (i + 1));
		}
		printf("All onEdges End\n");
		//return 0;
	}



	//printf("%s\n", newPathName);

	vector<int> nodes;
	nodes.reserve(g.n);
	for (int i = 0; i < g.n; i++) nodes.push_back(i);
	int* newNodes = new int[g.n];
	

	for (int i = 3; i < 4; i++)
	{
		int newN = 0.2 * (i + 1) * g.n;
		random_shuffle(nodes.begin(), nodes.end());

		for (int j = 0; j < newN; j++) newNodes[j] = nodes[j];

		Graph &subg = g.mksub(2, newNodes, newN);


		//char* cmdArr = new char[1000];
		FILE* fp;
		sprintf(newPathName, "%sNodes\\%s-0.%d-Node.txt", path, datasetsName, (i + 1) * 2);

		fp = fopen(newPathName, "w");

		fprintf(fp, "%d %d\n", subg.n, subg.e);

		for (int id = 0; id < newN; id++)
		{
			for (int j = subg.cd[id]; j < subg.cd[id] + subg.deg[id]; j++)
			{
				if(subg.adj[j] > id)
					fprintf(fp, "%d %d\n", id, subg.adj[j]);
			}
		}

		fclose(fp);

		printf("%lf Nodes End\n", 0.2 * (i + 1));
	}





	auto t3 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeGap(t1, t3)) / 1e6);

	printf("The End\n");

	return 0;

}
