// heap data structure :

template<class T>
class keyvalueLLU
{
public:
	int key;
	T value;
	keyvalueLLU& operator = (const keyvalueLLU& kv)
	{
		// if (*this == kv)
		//     return *this;
		this->key = kv.key;
		T tmp = kv.value;
		this->value = tmp;
		return *this;
	}
};

template<class T>
class bheapLLU
{
public:
	int n_max;// max number of nodes.
	int n;// number of nodes.
	int* pt;// pointers to nodes.
	keyvalueLLU<T>* kv;// (node,nck)
};

template<class T>
bheapLLU<T>* constructLLU(int n_max)
{
	bheapLLU<T>* heap = new bheapLLU<T>;
	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = new int[n_max];
	for (int i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = new keyvalueLLU<T>[n_max];
	return heap;
}
template<class T>
inline void swapLLU(bheapLLU<T>* heap, int i, int j)
{
	keyvalueLLU<T> kv_tmp = heap->kv[i];
	int pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}
template<class T>
inline void bubble_upLLU(bheapLLU<T>* heap, int i)
{
	int j = (i - 1) / 2;
	while (i > 0)
	{
		if (heap->kv[j].value > heap->kv[i].value)
		{
			swapLLU<T>(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

template<class T>
inline void bubble_downLLU(bheapLLU<T>* heap)
{
	int i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n)
	{
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value)
		{
			swapLLU<T>(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

template<class T>
inline void insertLLU(bheapLLU<T>* heap, keyvalueLLU<T> kv)
{
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_upLLU<T>(heap, heap->n - 1);
}

template<class T>
inline void updateLLU(bheapLLU<T>* heap, int key, T delta)
{
	int i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value) = delta;
		bubble_upLLU<T>(heap, i);
	}
}
template<class T>
inline keyvalueLLU<T> popminLLU(bheapLLU<T>* heap)
{
	keyvalueLLU<T> min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_downLLU(heap);
	return min;
}

template<class T>
inline keyvalueLLU<T> topLLU(bheapLLU<T>* heap)
{
	keyvalueLLU<T> min = heap->kv[0];
	return min;
}

template<class T>
//Building the heap structure with (key,value)=(node,k-clique degree) for each node
bheapLLU<T>* mkheapLLU(int n, T* nck)
{
	keyvalueLLU<T> kv;
	bheapLLU<T>* heap = constructLLU<T>(n);
	for (int i = 0; i < n; i++)
	{
		kv.key = i;
		kv.value = nck[i];
		insertLLU<T>(heap, kv);
	}
	return heap;
}
template<class T>
void freeheapLLU(bheapLLU<T>* heap)
{
	delete[] heap->pt;
	delete[] heap->kv;
	delete[] heap;
}