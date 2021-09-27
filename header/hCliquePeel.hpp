void hCliquePeeling(Graph& g, int h)
{
	int* GNodes = new int[g.n];
	for (int i = 0; i < g.n; i++) GNodes[i] = i;
	g.clique = new Clique(g.n, g.e, h);


	long long tol = 0;
	long long* cnt = new long long[g.n]();
	//g.kClique(h, &tol, cnt);

	g.kCliqueNew(h, &tol, cnt, GNodes, g.n);

	printf("Toltal cliques: %lld\n", tol);



	bheapLLU<long long>* cHeap = mkheapLLU<long long>(g.n, cnt);

	keyvalueLLU<long long> kv;

	long long curCliqueCore = 0, leftClique = tol;
	int leftN = g.n, leftM = g.e;
	int maxCliDenN = 0;	//the number of nodes in the subgraph achiving the largest hclique density;
	int maxCliDenM = 0;	//the number of edges in the subgraph achiving the largest hclique density;

	int maxCliCoreDenN = 0;	//the number of nodes in the maximal hclique core;
	int maxCliCoreDenM = 0;	//the number of edges in the maximal hclique core;
	double cliqueDensity = 0.0, curCliqueDensity, cliqueCoreDen = 0.0;

	int* nbrArr = new int[g.maxDeg], nbrNum;
	long long* nbrCnt = new long long[g.n](), nbrTol;
	long long maxCliDeg = 0;

	int coreTocore = 0, cntctc = 0;

	//peeling ordering infects maxCliqueDensity.
	while (leftN > 0)
	{

		kv = popminLLU<long long>(cHeap);
		long long cliqueDeg = kv.value;
		curCliqueDensity = 1.0 * leftClique / leftN;

		if (maxCliDeg < cliqueDeg)
		{
			cliqueCoreDen = curCliqueDensity;
			maxCliCoreDenN = leftN;
			maxCliCoreDenM = leftM;

			maxCliDeg = cliqueDeg;
			cntctc += coreTocore;
			//printf("maxCliDeg = %lld -> %lld, coreTocore = %d, \t leftN = %d, \t tolcn = %d\n", maxCliDeg - 1, maxCliDeg, coreTocore, leftN, cntctc);

			coreTocore = 0;
		}
		coreTocore++;

		maxCliDeg = max(maxCliDeg, cliqueDeg);

		//printf("id = %d, cliqueDeg = %lld, leftClique = %lld, leftN = %d\n", kv.key, kv.value, leftClique, leftN);



		//if (g.n - leftN < 100)
			//printf("-------------------------maxCliDeg = %lld, delId = %d, leftN = %d, leftClique: %lld\n", maxCliDeg, kv.key, leftN, leftClique);


		if (cliqueDensity < curCliqueDensity)
		{
			cliqueDensity = curCliqueDensity;
			maxCliDenN = leftN;
			maxCliDenM = leftM;
			//cliqueDensity = 1.0 * leftClique / leftN;
		}

		int delId = kv.key;
		nbrNum = 0;
		for (int i = g.cd[delId]; i < g.cd[delId] + g.deg[delId]; i++)
		{
			int adj = g.adj[i];
			nbrArr[nbrNum++] = adj;
		}
		nbrTol = 0;
		g.kCliqueNew(h - 1, &nbrTol, nbrCnt, nbrArr, nbrNum);


		leftClique -= nbrTol;



		for (int i = 0; i < nbrNum; i++)
		{
			int nbrId = nbrArr[i];
			cnt[nbrId] -= nbrCnt[nbrId];
			updateLLU(cHeap, nbrId, cnt[nbrId]);
			//printf("nbrId = %d, nbrCnt = %lld\n",nbrId, nbrCnt[nbrId]);
		}
		//printf("------------\n");


		for (int i = 0; i < nbrNum; i++)
			nbrCnt[nbrArr[i]] = 0;

		g.deleteNodes(&delId, 1);
		leftN--;
	}

	cntctc += coreTocore;
	printf("after coreTocore = %d, tolcn = %d\n", coreTocore, cntctc);

	printf("maxCliDenN = %d, maxCliDenM = %d, cliqueDensity = %lf\n", maxCliDenN, maxCliDenM, cliqueDensity);
	printf("maxCliDeg = %lld, maxCliCoreDenN = %d, cliqueCoreDensity = %lf\n", maxCliDeg, maxCliCoreDenN, cliqueCoreDen);
}