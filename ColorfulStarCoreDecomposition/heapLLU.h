// heap data structure :

class keyvalueLLU
{
public:
	int key;
	long long value;
};

class bheapLLU
{
public:
	int n_max;// max number of nodes.
	int n;// number of nodes.
	int* pt;// pointers to nodes.
	keyvalueLLU* kv;// (node,nck)
};

bheapLLU* constructLLU(int n_max)
{
	bheapLLU* heap = new bheapLLU;
	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = new int[n_max];
	for (int i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = new keyvalueLLU[n_max];
	return heap;
}

inline void swapLLU(bheapLLU * heap, int i, int j)
{
	keyvalueLLU kv_tmp = heap->kv[i];
	int pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

inline void bubble_upLLU(bheapLLU * heap, int i)
{
	int j = (i - 1) / 2;
	while (i > 0)
	{
		if (heap->kv[j].value > heap->kv[i].value)
		{
			swapLLU(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

inline void bubble_downLLU(bheapLLU * heap)
{
	int i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n)
	{
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value)
		{
			swapLLU(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

inline void insertLLU(bheapLLU * heap, keyvalueLLU kv)
{
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_upLLU(heap, heap->n - 1);
}

inline void updateLLU(bheapLLU * heap, int key, long long delta) 
{
	int i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value) = delta;
		bubble_upLLU(heap, i);
	}
}

inline keyvalueLLU popminLLU(bheapLLU * heap) 
{
	keyvalueLLU min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_downLLU(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,k-clique degree) for each node
bheapLLU* mkheapLLU(long long* nck, int n) 
{
	keyvalueLLU kv;
	bheapLLU* heap = constructLLU(n);
	for (int i = 0; i < n; i++) 
	{
		kv.key = i;
		kv.value = nck[i];
		insertLLU(heap, kv);
	}
	return heap;
}

void freeheapLLU(bheapLLU * heap)
{
	delete[] heap->pt;
	delete[] heap->kv;
	delete[] heap;
}