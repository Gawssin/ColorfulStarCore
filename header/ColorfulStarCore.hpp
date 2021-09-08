
void initColStarDegree(Graph& g, double** dp, int h, int colorNum, int* color, int** CC)
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

		//int* tol = new int[k + 1];

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
		}
		fill(NotColor0, NotColor0 + h, 0.0);
		fill(MustColor0, MustColor0 + h, 0.0);
	}
	delete[] NotColor0;
	delete[] MustColor0;
}

void ColorfulStarCoreDecomp(Graph& g, double** dp, int h, int* color, int** CC, double* ColofulStarCoreNum, double *maxCore = 0, int * maxCoreNum = 0)
{
	double* tmpDP = new double[g.n];
	for (int i = 0; i < g.n; i++) tmpDP[i] = dp[i][h - 1];
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
		kv = popminLLU<double>(heap);

		//printf("id = %d value = %lf\n", kv.key, kv.value);

		if (kv.value > starCoreNum)
		{
			starCoreNum = kv.value;
			maxN = leftN;
			maxM = leftM;
		}


		if (times % 50000 == 0)
			printf("times = %d left nodes = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", times, leftN, starCoreNum, maxN, maxM, 1.0 * maxM / maxN);

		//if (times == g.n - 2)
		//	printf("kv.value = %lf times = %d leftN = %d, leftM = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", kv.value, times, leftN, leftM, starCoreNum, maxN, maxM, 1.0 * maxM / maxN);


		//////
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
			updateLLU(heap, nbr, dp[nbr][h - 1]);
		}
		g.deg[i] = 0;
	}

	*maxCore = starCoreNum;
	*maxCoreNum = maxN;
	/////
	printf("End: times = %d left nodes = %d tolMax = %lf maxN = %d maxM = %d density = %lf\n", times, leftN, starCoreNum, maxN, maxM, 1.0 * maxM / maxN);

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
		//printf("%lf\n", tmpDP[i]);
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

		//delNodes[kv.key] = true;


		//printf("id = %d value = %lf\n", kv.key, kv.value);

		if (kv.value > starCoreNum)
		{
			starCoreNum = kv.value;
			maxN = leftN;
			maxM = leftM;
		}


		//////
		leftN--;
		int i = kv.key;
		leftM -= g.deg[i];
		//ColofulStarCoreNum[i] = starCoreNum;

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

	/////


	delete[] NotColor;
	delete[] MustColor;
}