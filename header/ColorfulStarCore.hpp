
void initColStarDegree(Graph& g, double** dp, int h, int colorNum, int* color, int** CC, int basicFlag)
{
	double* NotColor0 = new double[h]();
	double* MustColor0 = new double[h]();
	for (int i = 0; i < g.n; i++)
	{
		dp[i] = new double[h];
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

void ColorfulStarCoreDecomp(Graph& g, double** dp, int h, int* color, int** CC, double* ColofulStarCoreNum, int colorNum, int basicFlag, double *maxCore = 0, int * maxCoreNum = 0)
{
	double* tmpDP = new double[g.n];
	for (int i = 0; i < g.n; i++) tmpDP[i] = dp[i][h - 1];
	bheapLLU<double>* heap = mkheapLLU<double>(g.n, tmpDP);

	double maxStarDegree = -1;
	for (int i = 0; i < g.n; i++)
	{
		maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
	}
	printf("The maximum colorful h-star degree: %lf\n", maxStarDegree);



	int leftN = g.n, leftM = g.e;

	double* NotColor = new double[h]();
	double* MustColor = new double[h]();

	double starCoreNum = 0;
	int times = 0, maxN = 0, maxM = 0;
	keyvalueLLU<double> kv;

	printf("Times\t\tLeft_Nodes\tMaxCore(N)\tMaxCore(M)\tMaxCore(Density)\tMaxCore(CoreNumber)\n");

	while (leftN > 0)
	{
		times++;
		kv = popminLLU<double>(heap);

		if (kv.value > starCoreNum)
		{
			starCoreNum = kv.value;
			maxN = leftN;
			maxM = leftM;
		}

		//if (times % 100000 == 0)
			//printf("times = %d\tleft nodes = %-10d\ttolMax = %lf\tmaxN = %d\tmaxM = %d\tdensity = %lf\n", times, leftN, starCoreNum, maxN, maxM, maxN?(1.0 * maxM / maxN) : 0.0);

		if (times % 100000 == 0)
			printf("%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%lf\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, starCoreNum);

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

	if(maxCore != 0) *maxCore = starCoreNum;
	if(maxCoreNum != 0) *maxCoreNum = maxN;

	printf("End:\n%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%lf\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, starCoreNum);

	delete[] NotColor;
	delete[] MustColor;
}

void ColorfulStarCore(Graph& g, double** dp, int h, int* color, int** CC, double LB, int& delNum, int& delEdges)
{
	printf("LB: %lf\n", LB);

	double* tmpDP = new double[g.n];
	for (int i = 0; i < g.n; i++) tmpDP[i] = dp[i][h - 1];

	int sdp = 0;
	for (int i = 0; i < g.n; i++)
	{
		if (tmpDP[i] == 0.0) sdp++;
	}

	printf("sdp = %d\n", sdp);

	bheapLLU<double>* heap = mkheapLLU<double>(g.n, tmpDP);

	double maxStarDegree = -1;
	for (int i = 0; i < g.n; i++)
	{
		maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
	}
	printf("maxStarDegree = %lf\n", maxStarDegree);

	int leftN = g.n, leftM = g.e;

	double* NotColor = new double[h]();
	double* MustColor = new double[h]();

	double starCoreNum = 0;
	int times = 0, maxN = 0, maxM = 0;
	keyvalueLLU<double> kv;

	while (leftN > 0)
	{
		times++;
		double Min = 1e300;
		kv = popminLLU<double>(heap);

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
			updateLLU<double>(heap, nbr, dp[nbr][h - 1]);
		}
		g.deg[i] = 0;

	}

	delete[] NotColor;
	delete[] MustColor;
}