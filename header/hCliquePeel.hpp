void hCliquePeeling(Graph& g, int h)
{
	int* GNodes = new int[g.n];
	for (int i = 0; i < g.n; i++) GNodes[i] = i;
	g.clique = new Clique(g.n, g.e, h);

	long long tol = 0;
	long long* cnt = new long long[g.n]();

	g.kCliqueNew(h, &tol, cnt, GNodes, g.n);

	printf("Total cliques: %lld\n", tol);

	bheapLLU<long long>* cHeap = mkheapLLU<long long>(g.n, cnt);

	keyvalueLLU<long long> kv;

	long long curCliqueCore = 0, leftClique = tol;
	int leftN = g.n, leftM = 0;

	for (int i = 0; i < g.n; i++)
		leftM += g.deg[i];
	leftM /= 2;

	int maxCliDenN = 0;	//the number of nodes in the subgraph achiving the largest h-clique density;
	int maxCliDenM = 0;	//the number of edges in the subgraph achiving the largest h-clique density;

	int maxCliCoreDenN = 0;	//the number of nodes in the h-clique Kmax core;
	int maxCliCoreDenM = 0;	//the number of edges in the h-clique Kmax core;
	double cliqueDensity = 0.0, curCliqueDensity, cliqueCoreDen = 0.0;

	int* nbrArr = new int[g.maxDeg], nbrNum;
	long long* nbrCnt = new long long[g.n](), nbrTol;
	long long maxCliDeg = 0;

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
		}

		maxCliDeg = max(maxCliDeg, cliqueDeg);

		if (cliqueDensity < curCliqueDensity)
		{
			cliqueDensity = curCliqueDensity;
			maxCliDenN = leftN;
			maxCliDenM = leftM;
		}

		int delId = kv.key;
		leftM -= g.deg[delId];
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
		}

		for (int i = 0; i < nbrNum; i++)
			nbrCnt[nbrArr[i]] = 0;

		g.deleteNodes(&delId, 1);
		leftN--;
	}

	printf("\nH-clique densest subgraph (approximation)\n");
	printf("Nodes:\t\t%d\nEdges:\t\t%d\nClique-Density:\t%lf\n\n\n", maxCliDenN, maxCliDenM, cliqueDensity);

	printf("H-clique Kmax core\n");
	printf("Nodes:\t\t%d\nEdges:\t\t%d\nClique-Density:\t%lf\nKmax:\t\t%d\n\n", maxCliCoreDenN, maxCliCoreDenM, cliqueCoreDen, maxCliDeg);
}