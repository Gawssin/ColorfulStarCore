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
using namespace std;

#define debug 1

long long kCCN = 0, *tmpCnt, *cntSub;
int* uadj, * uadjt, * Merge, * mark, * oloc;


int main(int argc, char** argv)
{
	auto t0 = getTime();

	Graph g;
	int k = atoi(argv[1]);
	cout << "Reading edgelist from file " << argv[2] << endl;
	g.readedgelist(argv[2]);
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

	//double*** ck = new double** [g.n];

	printf("colorNum = %d\n", colorNum);

	int** CC = new int* [g.n];

	//int** nbrCol = new int* [g.n];


	for (int i = 0; i < g.n; i++)
	{
		dp[i] = new double[k];
		int colorNum_i = 0;
		//int* C = new int[colorNum]();
		CC[i] = new int[colorNum]();
		for (int j = g.cd[i]; j < g.cd[i + 1]; j++)
		{
			int nbr = g.adj[j];
			CC[i][color[nbr]]++;
			colorNum_i = max(colorNum_i, color[nbr]);
		}
		double* NotColor0 = new double[k]();
		double* MustColor0 = new double[k]();
		//int* tol = new int[k + 1];

		NotColor0[0] = 1;
		for (int c = 1; c <= colorNum_i; c++)
		{
			for (int j = k - 1; j > 0; j--)
				NotColor0[j] = NotColor0[j - 1] * CC[i][c] + NotColor0[j];
		}

		for (int j = 1; j < k; j++)
		{
			MustColor0[j] = NotColor0[j - 1] * CC[i][0];
			dp[i][j] = MustColor0[j] + NotColor0[j];
		}
		delete[] NotColor0;
		delete[] MustColor0;

		//------------------------------------------------------------------

		//delete[] C;

		//if (i == 2660)
		//{
		//	printf("%d %d \n", cNum, colorNum_i);
		//	for (int q = 0; q < colorNum_i; q++)
		//	{
		//		printf("%d ", C[q]);

		//	}
		//	printf("\n");
		//	for (int q = 0; q < colorNum_i; q++)
		//	{
		//		printf("%d ", nbrCol[i][q]);

		//	}
		//	printf("\n");
		//if (dp[i][k] < 0)



		//printf("%d %lf %lf\n", i, dp[i][k - 1], ck[cNum][k - 1]);

		//if (dp[i][k] != ck[i][cNum][k])
			//printf("Not %d colorNum = %d %lf %lf\n", i, colorNum_i, dp[i][k], ck[i][cNum][k]);

		//}
	}

	auto t2 = getTime();

	bheapLLU* heap = mkheapLLU(g.n, dp, k - 1);

	if (debug)
	{
		for (int i = 0; i < g.n; i++)
		{
			printf("id = %d starDegree = %lf\n", i, dp[i][k - 1]);
		}
	}

	double maxStarDegree = -1;
	for (int i = 0; i < g.n; i++)
	{
		maxStarDegree = max(maxStarDegree, dp[i][k - 1]);
		//printf("id = %d starDegree = %lf\n", i, dp[i][k - 1]);
	}


	printf("maxStarDegree = %lf\n", maxStarDegree);


	int delNodes = 0, times = 0;
	bool* mark = new bool[g.n]();
	int leftN = g.n, leftM = g.e;

	double* NotColor0 = new double[k]();
	double* MustColor0 = new double[k]();
	double* ColofulStarCoreNum = new double[g.n];
	//int* C = new int[colorNum]();

	double tolMax = 0;
	int maxN = 0, maxM = 0;

	keyvalueLLU kv;
	while (leftN > 0)
	{
		times++;
		double Min = 1e300;
		kv = popminLLU(heap);

		int revId = kv.key;
		//printf("id = %d value = %lf\n", kv.key, kv.value);
		Min = min(Min, kv.value);

		if (Min > tolMax)
		{
			tolMax = Min;
			maxN = leftN;
			maxM = leftM;
			//printf("tolMax = %lf maxN = %d\n", tolMax, maxN);
		}
		//tolMax = max(tolMax, Min);

		if (times % 50000 == 0)
			printf("times = %d left nodes = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", times, leftN, tolMax, maxN, maxM, 1.0 * maxM / maxN);

		leftN--;

		//////

		int i = revId;
		leftM -= g.deg[i];
		ColofulStarCoreNum[i] = tolMax;

		for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
		{
			int nbr = g.adj[j];

			for (int h = g.cd[nbr]; h < g.cd[nbr] + g.deg[nbr]; h++)
			{
				int hnbr = g.adj[h];
				if (hnbr == i)
				{
					swap(g.adj[h], g.adj[g.cd[nbr] + g.deg[nbr] - 1]);
					g.deg[nbr]--;
					break;
				}
			}


			CC[nbr][color[i]]--;


			double lastdp = dp[nbr][k - 1];
			NotColor0[0] = 1;
			for (int h = 1; h < k; h++)
			{
				MustColor0[h] = NotColor0[h - 1] * (CC[nbr][color[i]] + 1);
				NotColor0[h] = dp[nbr][h] - MustColor0[h];
			}

			for (int h = 1; h < k; h++)
			{
				MustColor0[h] = NotColor0[h - 1] * CC[nbr][color[i]];
				dp[nbr][h] = NotColor0[h] + MustColor0[h];
			}
			updateLLU(heap, nbr, dp[nbr][k - 1]);

		}

	}

	/////

	printf("End: times = %d left nodes = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", times, leftN, tolMax, maxN, maxM, 1.0 * maxM / maxN);

	if (debug)
	{
		for (int i = 0; i < g.n; i++)
		{
			printf("id = %d starDegree = %lf\n", i, ColofulStarCoreNum[i]);
		}
	}


	delete[] NotColor0;
	delete[] MustColor0;
	//delete[] C;


	auto t3 = getTime();
	printf("- Overall time = %lfs\n", ((double)timeDrt(t2, t3)) / 1e6);

	printf("The End\n");

	return 0;

}
