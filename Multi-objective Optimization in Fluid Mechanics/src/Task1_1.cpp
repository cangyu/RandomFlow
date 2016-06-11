/**
 * Task1_1,find the minimum value of the given single variable function using genetic algorithm.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cstring>

using namespace std;

const double eps = 1e-12;

double sectionLeft = -5.0;
double sectionRight = 5.0;

double precision = 1e-5;
size_t PopulationSize = 200;
size_t GenerationCnt = 40;
double P_Cross = 0.65;
double P_Mutate = 0.001;

double objFunc(double x)
{
	return pow(x - 2.5, 2);
}

class chromosome
{
private:
	string gene;
	double realVal, objVal, fitVal;

public:
	chromosome(){}
	chromosome(size_t len) :gene(len, '0')
	{
		randomize();
		decode();
		calcObj();
		calcFit();
	}
	chromosome(const chromosome &rhs) :gene(rhs.gene), realVal(rhs.realVal), objVal(rhs.objVal), fitVal(rhs.fitVal){}
	chromosome& operator=(chromosome rhs)
	{
		swap(gene, rhs.gene);
		swap(realVal, rhs.realVal);
		swap(objVal, rhs.objVal);
		swap(fitVal, rhs.fitVal);
		return *this;
	}
	~chromosome() {}

	bool operator<(const chromosome &rhs) const { return fitVal < rhs.fitVal; }

	/*产生随机的基因序列*/
	void randomize()
	{
		for (size_t i = 0; i < gene.length(); i++)
		{
			int choice = rand() % 1000;
			if (choice < 500)
				gene[i] = '0';
			else
				gene[i] = '1';
		}
	}

	/*将二进制基因序列转换成十进制的小数*/
	double decode()
	{
		double tmp = 0;
		const int n = gene.length();
		for (int i = n - 1; i >= 0; i--)
			if (gene[i] == '1')
				tmp += pow(2, n - 1 - i);

		return realVal = min(sectionLeft + tmp*precision, sectionRight);
	}

	/*计算目标函数值*/
	double calcObj()
	{
		return objVal = objFunc(realVal);
	}

	/*计算适应值*/
	double calcFit()
	{
		return fitVal = -objVal;
	}

	/*重新计算所有值*/
	void update()
	{
		decode();
		calcObj();
		calcFit();
	}

	/*获取相关结果*/
	double getRealVal() const { return realVal; }
	double getObjVal() const { return objVal; }
	double getFitVal() const { return fitVal; }

	/*随机单点交叉*/
	void crossOver(chromosome &partner, double probability)
	{
		if (gene.length() != partner.gene.length())
			throw "Chromosome not match!\n";
		if (probability < 0 || probability>1)
			throw "Invalid probability!\n";

		bool doSwap = false;
		if (fabs(probability) < eps)//预设的交叉概率为0，保持不变
			return;
		else if (fabs(probability - 1) < eps)//预设的变异概率为1，必然有部分要进行交换
			doSwap = true;
		else
		{
			double choice = (rand() % 32749) / 32749.0;
			doSwap |= choice < probability;
		}

		if (doSwap)
		{
			int pos = rand() % gene.length();//Position to start swapping
			for (int i = pos; i < gene.length(); i++)
				swap(gene[i], partner.gene[i]);
		}
	}

	/*二进制位变异*/
	void mutate(double probability)
	{
		if (probability < 0 || probability>1)
			throw "Invalid probability!\n";

		const int total = '0' + '1';

		if (fabs(probability) < eps)//预设的变异概率为0，保持不变
			return;
		else if (fabs(probability - 1) < eps)//预设的变异概率为1，全部取反
			for (int i = 0; i < gene.length(); i++)
				gene[i] = total - gene[i];
		else
		{
			double choice = 0;
			for (int i = 0; i < gene.length(); i++)
			{
				choice = (rand() % 32749) / 32749.0;
				if (choice < probability)
					gene[i] = total - gene[i];
			}
		}
	}
};

int main(int argc, char **argv)
{
	ofstream fout("./results/Test1_1.txt");
	if (!fout)
		throw "Invalid output file path!\n";

	//输入参数
	cout << "目标函数：f(x)=(x-2.5)^2" << endl;
	cout << "请输入相关参数：" << endl;
	cout << "解区间(以空格间开，e.g. -5 5)："; cin >> sectionLeft >> sectionRight;
	cout << "解的精度(e.g. 0.001 or 1e-3)："; cin >> precision;
	cout << "种群规模(e.g. 200 or 2e2)："; cin >> PopulationSize;
	cout << "迭代次数(e.g. 50 or 5e1)："; cin >> GenerationCnt;
	cout << "交叉概率(e.g. 0.65 or 6.5e-1)："; cin >> P_Cross;
	cout << "变异概率(e.g. 0.001 or 1e-3)："; cin >> P_Mutate;
	cout << "Calculating..." << endl;

	//圆整
	size_t chromoLen = ceil(log2f((sectionRight - sectionLeft) / precision));
	precision = (sectionRight - sectionLeft) / pow(2, chromoLen);

	fout << "当前测试参数如下：" << endl;
	fout << "目标函数：f(x)=(x-2.5)^2" << endl;
	fout << "取值区间：[" << sectionLeft << "," << sectionRight << "]" << endl;
	fout << "圆整后的实际精度：" << precision << endl;
	fout << "基因串长度：" << chromoLen << endl;
	fout << "种群规模：" << PopulationSize << endl;
	fout << "迭代次数：" << GenerationCnt << endl;
	fout << "交叉概率：" << P_Cross << endl;
	fout << "变异概率：" << P_Mutate << endl;
	fout << endl;

	//初始化种群
	srand(time(NULL));
	vector<chromosome> grp_cur;
	for (size_t i = 0; i < PopulationSize; i++)
		grp_cur.push_back(chromosome(chromoLen));
	sort(grp_cur.begin(), grp_cur.end());

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

			grp_next[c1].crossOver(grp_next[c2], P_Cross);
			hasCrossed[c1] = true;
			hasCrossed[c2] = true;
		}

		//变异
		for (size_t i = 0; i < PopulationSize; i++)
			grp_next[i].mutate(P_Mutate);

		//Update
		for (size_t i = 0; i < PopulationSize; i++)
			grp_next[i].update();
		sort(grp_next.begin(), grp_next.end());

		//Keep the best
		//swap(grp_next[0], grp_cur[grp_cur.size() - 1]);
		//sort(grp_next.begin(), grp_next.end());
		grp_cur = grp_next;

		fout << "第" << k << "次迭代：" << endl
			 << "最优值为：" << grp_cur.back().getRealVal()
			 << "，在目标函数作用下的结果为：" << grp_cur.back().getObjVal()
			 << "，适应度为：" << grp_cur.back().getFitVal() << endl;
			 //<< "最差值为：" << grp_cur.front().getRealVal()
			 //<< " ，在目标函数作用下的结果为：" << grp_cur.front().getObjVal()
			 //<< " ,适应度为：" << grp_cur.front().getFitVal() << endl;
	}

	delete[] hasCrossed;
	fout.close();
	return 0;
}

