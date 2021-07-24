#include <chrono>
using namespace chrono;

long long combination(int n, int m)
{
	if (n < m) return 0;
	if (n == m) return 1;
	long long res = 1;
	for (long long i = n, j = 1; i >= n - m + 1; i--, j++) res *= i, res /= j;
	return res;
}

int outList1(int* list1, int* i, int* j, int s, int s1)
{
	if (*i == s - 1)
		return list1[++ * j];
	if (*j == s1 - 1)
		return list1[++ * i];
	if (list1[*i + 1] < list1[*j + 1])
		return list1[++ * i];
	return list1[++ * j];
}

int merging(int s, int* list1, int s1, int* list2, int s2, int* list3)
{
	int i = 0, j = 0, p = -1, q = s - 1, s3 = 0;
	int x = outList1(list1, &p, &q, s, s1), y = list2[0];
	while (i < s1 && j < s2)
	{
		if (x < y)
		{
			x = outList1(list1, &p, &q, s, s1);
			++i;
			//x = list1[++i];
			continue;
		}
		if (y < x)
		{
			y = list2[++j];
			continue;
		}
		list3[s3++] = x;
		x = outList1(list1, &p, &q, s, s1);
		++i;
		//x = list1[++i];
		y = list2[++j];
	}
	return s3;
}


long long Binary(double left, double right, Graph& g, double* uBd)
{
	//int left = 0, right = n;      //解的范围初始为(0,n],不包含0                             
	while (left + 1 < right)
	{
		long long mid = (right + left) / 2;
		if (1)//check()
		{
			right = mid;          //修正解的范围(left,right]
		}
		else
			left = mid;
	}
	return right;                //最后left + 1 = right
}

high_resolution_clock::time_point getTime()
{
	return high_resolution_clock::now();
}


auto timeGap(high_resolution_clock::time_point t1, high_resolution_clock::time_point t2)
{
	return duration_cast<microseconds>(t2 - t1).count();
}