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


high_resolution_clock::time_point getTime()
{
	return high_resolution_clock::now();
}


auto timeGap(high_resolution_clock::time_point t1, high_resolution_clock::time_point t2)
{
	return duration_cast<microseconds>(t2 - t1).count();
}




class readCMD
{
	int argc;
	char** argv;
	int index;
	char* arr;
public:
	readCMD(int argc, char** argv) : argc(argc), argv(argv), index(0) { arr = new char[1000]; };
	~readCMD() { delete[] arr; };
	char* read();
};

char* readCMD::read()
{
	index++;
	char* cmdArr = new char[1000];
	FILE* fp;
	fp = fopen(argv[1], "r");
	for (int i = 0; i < index; i++)
	{
		fscanf(fp, "%s", cmdArr);
	}
	fclose(fp);
	return cmdArr;
}

void _int128_print(__int128 x)
{
	if (!x) return;
	if (x < 0) putchar('-'), x = -x;
	_int128_print(x / 10);
	putchar(x % 10 + '0');
}


char* _int128_to_str(__int128 x)
{
	if (x / 10 == 0)
	{
		char* str = new char[40]();
		if (x < 0)
		{
			str[0] = '-';
			str[1] = -x + '0';
		}
		else
		{
			str[0] = x + '0';
		}
		return str;
	}

	char* str = _int128_to_str(x / 10);

	if (x < 0) x = -x;
	str[strlen(str)] = x % 10 + '0';
	return str;
}