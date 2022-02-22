
void initColStarDegree(Graph& g, __int128** dp, int h, int colorNum, int* color, int** CC, int basicFlag)
{
	__int128* NotColor0 = new __int128[h]();
	__int128* MustColor0 = new __int128[h]();
	for (int i = 0; i < g.n; i++)
	{
		dp[i] = new __int128[h];
		int colorNum_i = 0;
		//int* C = new int[colorNum]();
		CC[i] = new int[colorNum]();
		for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
		{
			int nbr = g.adj[j];
			CC[i][color[nbr]]++;
			colorNum_i = max(colorNum_i, color[nbr]);
		}

		if (basicFlag == 0)
		{
			NotColor0[0] = 1;
			for (int c = 1; c <= colorNum_i; c++)
			{
				for (int j = h - 1; j > 0; j--)
					NotColor0[j] = NotColor0[j - 1] * CC[i][c] + NotColor0[j];
			}

			for (int j = 1; j < h; j++)
			{
				MustColor0[j] = NotColor0[j - 1] * CC[i][0];
				dp[i][j] = MustColor0[j] + NotColor0[j];
				NotColor0[j - 1] = 0;
			}
			NotColor0[h - 1] = 0;
		}
		else
		{
			dp[i][0] = 1;
			for (int c = 0; c <= colorNum_i; c++)
			{
				for (int j = h - 1; j > 0; j--)
					dp[i][j] = dp[i][j] + dp[i][j - 1] * CC[i][c];
			}
		}
	}
	delete[] NotColor0;
	delete[] MustColor0;
}

void ColorfulStarCoreDecomp(Graph& g, __int128** dp, int h, int* color, int** CC, __int128* ColofulStarCoreNum, int colorNum, int basicFlag, __int128* maxCore = 0, int* maxCoreNum = 0)
{
	__int128* tmpDP = new __int128[g.n];
	for (int i = 0; i < g.n; i++) tmpDP[i] = dp[i][h - 1];
	bheapLLU<__int128>* heap = mkheapLLU<__int128>(g.n, tmpDP);

	__int128 maxStarDegree = -1;
	for (int i = 0; i < g.n; i++)
	{
		maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
	}
	printf("The maximum colorful h-star degree: %s\n", _int128_to_str(maxStarDegree));



	int leftN = g.n, leftM = g.e;

	__int128* NotColor = new __int128[h]();
	__int128* MustColor = new __int128[h]();

	__int128 starCoreNum = 0;
	int times = 0, maxN = 0, maxM = 0;
	keyvalueLLU<__int128> kv;

	printf("Times\t\tLeft_Nodes\tMaxCore(N)\tMaxCore(M)\tMaxCore(Density)\tMaxCore(CoreNumber)\n");

	while (leftN > 0)
	{
		times++;
		kv = popminLLU<__int128>(heap);

		if (kv.value > starCoreNum)
		{
			starCoreNum = kv.value;
			maxN = leftN;
			maxM = leftM;
		}

		//if (times % 100000 == 0)
			//printf("times = %d\tleft nodes = %-10d\ttolMax = %lf\tmaxN = %d\tmaxM = %d\tdensity = %lf\n", times, leftN, starCoreNum, maxN, maxM, maxN?(1.0 * maxM / maxN) : 0.0);

		if (times % 100000 == 0)
			printf("%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%s\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, _int128_to_str(starCoreNum));

		leftN--;
		int i = kv.key;
		leftM -= g.deg[i];
		ColofulStarCoreNum[i] = starCoreNum;

		for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
		{
			int nbr = g.adj[j];

			for (int p = g.cd[nbr]; p < g.cd[nbr] + g.deg[nbr]; p++)
			{
				int hnbr = g.adj[p];
				if (hnbr == i)
				{
					swap(g.adj[p], g.adj[g.cd[nbr] + g.deg[nbr] - 1]);
					g.deg[nbr]--;
					break;
				}
			}

			CC[nbr][color[i]]--;

			NotColor[0] = 1;
			if (basicFlag == 0)
			{
				for (int p = 1; p < h; p++)
				{
					MustColor[p] = NotColor[p - 1] * (CC[nbr][color[i]] + 1);
					NotColor[p] = dp[nbr][p] - MustColor[p];
				}
				for (int p = 1; p < h; p++)
				{
					MustColor[p] = NotColor[p - 1] * CC[nbr][color[i]];
					dp[nbr][p] = NotColor[p] + MustColor[p];
				}
			}
			else
			{
				fill(dp[nbr] + 1, dp[nbr] + h, 0);
				for (int c = 0; c < colorNum; c++)
				{
					for (int j = h - 1; j > 0; j--)
						dp[nbr][j] = dp[nbr][j] + dp[nbr][j - 1] * CC[nbr][c];
				}
			}

			updateLLU(heap, nbr, dp[nbr][h - 1]);
		}
		g.deg[i] = 0;
	}

	if (maxCore != 0)* maxCore = starCoreNum;
	if (maxCoreNum != 0)* maxCoreNum = maxN;

	printf("End:\n%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%s\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, _int128_to_str(starCoreNum));

	delete[] NotColor;
	delete[] MustColor;
}

void ColorfulStarCore(Graph& g, __int128** dp, int h, int* color, int** CC, __int128 LB, int& delNum, int& delEdges)
{
	printf("Lower Bound: %s\n", _int128_to_str(LB));

	__int128* tmpDP = new __int128[g.n];
	for (int i = 0; i < g.n; i++) tmpDP[i] = dp[i][h - 1];

	bheapLLU<__int128>* heap = mkheapLLU<__int128>(g.n, tmpDP);

	__int128 maxStarDegree = -1;
	for (int i = 0; i < g.n; i++)
	{
		maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
	}
	printf("maxStarDegree = %s\n", _int128_to_str(maxStarDegree));

	int leftN = g.n, leftM = g.e;

	__int128* NotColor = new __int128[h]();
	__int128* MustColor = new __int128[h]();

	__int128 starCoreNum = 0;
	int times = 0, maxN = 0, maxM = 0;
	keyvalueLLU<__int128> kv;

	while (leftN > 0)
	{
		times++;
		//__int128 Min = 1e300;
		kv = popminLLU<__int128>(heap);

		if (kv.value >= LB)
		{
			break;
		}
		delNum++;
		delEdges += g.deg[kv.key];

		if (kv.value > starCoreNum)
		{
			starCoreNum = kv.value;
			maxN = leftN;
			maxM = leftM;
		}

		leftN--;
		int i = kv.key;
		leftM -= g.deg[i];

		for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
		{
			int nbr = g.adj[j];

			for (int p = g.cd[nbr]; p < g.cd[nbr] + g.deg[nbr]; p++)
			{
				int hnbr = g.adj[p];
				if (hnbr == i)
				{
					swap(g.adj[p], g.adj[g.cd[nbr] + g.deg[nbr] - 1]);
					g.deg[nbr]--;
					break;
				}
			}

			CC[nbr][color[i]]--;

			NotColor[0] = 1;
			for (int p = 1; p < h; p++)
			{
				MustColor[p] = NotColor[p - 1] * (CC[nbr][color[i]] + 1);
				NotColor[p] = dp[nbr][p] - MustColor[p];
			}

			for (int p = 1; p < h; p++)
			{
				MustColor[p] = NotColor[p - 1] * CC[nbr][color[i]];
				dp[nbr][p] = NotColor[p] + MustColor[p];
			}
			updateLLU<__int128>(heap, nbr, dp[nbr][h - 1]);
		}
		g.deg[i] = 0;

	}

	delete[] NotColor;
	delete[] MustColor;
}