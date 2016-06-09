/**
 * Task1_1,find the minimum value of the given single variable function using genetic algorithm.
 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

double sectionLeft = -5.0;
double sectionRight = 5.0;

double precision = 1e-6;
size_t PopulationSize = 1000;
size_t GenerationCnt = 200;
double P_Cross = 0.7;
double P_Mutate = 0.01;

bool checkParemeters();
double objFunc(double x);

class chromosome
{
private:
	string gene;
	double objVal, fitVal;

public:
	chromosome(){}
	chromosome(size_t len) :gene(len, '0'), objVal(0), fitVal(0)
	{
		randomize();
	}
	chromosome(const chromosome &rhs) :gene(rhs.gene), objVal(rhs.objVal), fitVal(rhs.fitVal){}
	chromosome& operator=(chromosome rhs)
	{
		swap(gene, rhs.gene);
		swap(objVal, rhs.objVal);
		swap(fitVal, rhs.fitVal);
		return *this;
	}
	~chromosome() {}

	bool operator<(const chromosome &rhs){ return fitVal < rhs.fitVal; }

	/*产生随机的基因序列*/
	void randomize()
	{}

	/*将基因串清0*/
	void reset()
	{}

	double getFitVal() const { return fitVal; }

	/*随机单点交叉*/
	void cross_over(chromosome &partner)
	{
		if (gene.length() != partner.gene.length())
			throw "Chromosome not match!\n";
	}


	/*二进制位变异*/
	void mutate()
	{

	}

	/*更新适应值*/
	int calcFit()
	{

		return fitVal;
	}


};

int main(int argc, char **argv)
{
	//输入参数
	do{
		cout << "请输入相关参数：" << endl;
		cout << "解区间(以空格间开，e.g. -5 5)："; cin >> sectionLeft >> sectionRight;
		cout << "解的精度(e.g. 0.001 or 1e-3)："; cin >> precision;
		cout << "种群规模(e.g. 200 or 2e2)："; cin >> PopulationSize;
		cout << "迭代次数(e.g. 50 or 5e1)："; cin >> GenerationCnt;
		cout << "交叉概率(e.g. 0.65 or 6.5e-1)："; cin >> P_Cross;
		cout << "变异概率(e.g. 0.001 or 1e-3)："; cin >> P_Mutate;
	} while (!checkParemeters());

	//圆整
	size_t chromoLen = ceil(log2f((sectionRight - sectionLeft) / precision));
	precision = (sectionRight - sectionLeft) / pow(2, chromoLen);
	cout << "圆整后的实际精度为：" << precision << endl;

	//初始化种群
	srand(time(NULL));
	vector<chromosome> grp_cur(PopulationSize, chromosome(chromoLen));

	//迭代进化
	bool *hasCrossed = new bool[PopulationSize];

	for (size_t k = 0; k < GenerationCnt; k++)
	{
		vector<chromosome> grp_next;		
		int c1 = -1, c2 = -1;

		//选择
		for (size_t cnt = 0; cnt < PopulationSize; cnt++)
		{
			do{
				c1 = rand() % PopulationSize;
				c2 = rand() % PopulationSize;
			} while (c1 == c2);
			
			double f1 = grp_cur[c1].getFitVal();
			double f2 = grp_cur[c2].getFitVal();

			if (f1 > f2)
				grp_next.push_back(grp_cur[c1]);
			else
				grp_next.push_back(grp_cur[c2]);
		}

		//交叉
		memset(hasCrossed, 0, PopulationSize*sizeof(bool));

		for (size_t i = 0; i < PopulationSize; i += 2)
		{
			do{
				c1 = rand() % PopulationSize;
				c2 = rand() % PopulationSize;
			} while (c1 == c2 || hasCrossed[c1] || hasCrossed[c2]);

			grp_next[c1].cross_over(grp_next[c2]);
			hasCrossed[c1] = true;
			hasCrossed[c2] = true;
		}

		//变异
		for (auto t : grp_next)
			t.mutate();

		//Update
		for (auto t : grp_next)
			t.calcFit();
		sort(grp_next.begin(), grp_next.end());

		//Keep the best


	}




	return 0;
}

bool checkParemeters()
{
	return true;
}

double objFunc(double x) { return pow(x - 2.5, 2); }
