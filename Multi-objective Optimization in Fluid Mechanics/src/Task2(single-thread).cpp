/**
 * Task2:Nash GA for minimization of F(x,y) with 2 players.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <limits>

using namespace std;

typedef double(*ObjFuncPtr)(const vector<double> &);

ofstream fout("../results/Test2(single-thread).txt");
const double eps = 1e-12;

size_t varCnt = 2;
vector<double> leftBound, rightBound, resolution;
vector<size_t> varBitLen;

size_t RoundCnt = 100;
vector<double> ans(2, 0);

size_t PopulationSize = 200;
size_t GenerationCnt = 40;
double P_Cross = 0.65;
double P_Mutate = 0.001;

/**
 * Calculate the object function of F(x,y)=((x-1)y)^2+(xy-2)^4
 */
double F(const vector<double> &var)
{
	return pow(var[0] - 1, 2)*pow(var[1], 2) + pow(var[0] * var[1] - 2, 4);
}

/**
 * Class definition for general chromosome used in this algorithm
 */
class chromosome
{
private:
	string gene;
	size_t targetIndex;
	ObjFuncPtr f;
	double objVal, fitVal;

	//swap the data members between 2 chromosomes
	void exchange(chromosome &rhs)
	{
		swap(gene, rhs.gene);
		swap(targetIndex, rhs.targetIndex);
		swap(f, rhs.f);
		swap(objVal, rhs.objVal);
		swap(fitVal, rhs.fitVal);
		swap(var, rhs.var);
	}

	//view the gene as encoded in gray code
	//xfer the gray code to its binary representation
	string gray2binary(const string &gstr)
	{
		string tmp(gstr.length(), '0');
		tmp[0] = gstr[0];

		for (size_t i = 1; i < tmp.length(); i++)
			if (tmp[i - 1] != gstr[i])
				tmp[i] = '1';

		return tmp;
	}

public:
	vector<double> var;

	//Constructors
	chromosome(){}
	chromosome(const ObjFuncPtr &_f, const vector<double> &_var, size_t _self, size_t _len) :gene(_len, '0'), targetIndex(_self), f(_f), var(_var), objVal(0), fitVal(0)
	{
		randomize();
		decode();
		calcObj();
		calcFit();
	}
	chromosome(const chromosome &rhs) :gene(rhs.gene), targetIndex(rhs.targetIndex), var(rhs.var), f(rhs.f), objVal(rhs.objVal), fitVal(rhs.fitVal){}
	
	//Destructor
	~chromosome() {}
	
	//operator=,using the copy-and-swap strategy
	chromosome& operator=(chromosome rhs)
	{
		exchange(rhs);
		return *this;
	}

	//operator<, which is used to sort the chromosomes in a group
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

	/*基因的二进制表示转换为十进制浮点数*/
	void decode()
	{
		string bstr = gray2binary(gene);
		double tmp = 0;
		const size_t n = bstr.length();
		for (size_t j = 0; j < n; j++)
			if (bstr[j] == '1')
				tmp += pow(2, n - 1 - j);
		var[targetIndex] = min(leftBound[targetIndex] + tmp*resolution[targetIndex], rightBound[targetIndex]);
	}

	/*计算目标函数值*/
	double calcObj() { return objVal = (*f)(var); }

	/*计算适应值*/
	double calcFit() { return fitVal = -objVal; }

	/*重新计算所有值*/
	void update()
	{
		decode();
		calcObj();
		calcFit();
	}

	/*获取相关结果*/
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
			double choice = (rand() % 100) / 100.0;
			doSwap |= choice < probability;
		}

		if (doSwap)
		{
			int pos = rand() % gene.length();//Position to start swapping
			for (size_t i = pos; i < gene.length(); i++)
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
			for (size_t i = 0; i < gene.length(); i++)
				gene[i] = total - gene[i];
		else
		{
			double choice = 0;
			for (size_t i = 0; i < gene.length(); i++)
			{
				choice = (rand() % 100) / 100.0;
				if (choice < probability)
					gene[i] = total - gene[i];
			}
		}
	}
};

/**
 * Utilizing the classic genetic algorithm to optimize a given variable in the object function.
 *
 * The selection operator is implemented by the tournament strategy;
 * The cross-over operator is implemented by crossing over a random index with a given possibility P_Cross to decide whether the operation will be taken;
 * The mutation operator is implemented as usual with a given possibility P_Mutate to decide whether the operation will be taken on each bit.
 */
void GA_Solver(ObjFuncPtr func, vector<double> &var, size_t self)
{
	//初始化种群
	vector<chromosome> grp_cur;
	for (size_t i = 0; i < PopulationSize; i++)
		grp_cur.push_back(chromosome(func, var, self, varBitLen[self]));
	sort(grp_cur.begin(), grp_cur.end());

	//迭代进化
	vector<bool> hasCrossed(PopulationSize, false);

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
		fill(hasCrossed.begin(), hasCrossed.end(), false);

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
	}

	//返回最优
	var[self] = grp_cur.back().var[self];
}

/**
 * Main body that performs the Nash progress.
 */
int main(int argc, char **argv)
{
	if (!fout)
		throw "Invalid output file path!\n";
	
	srand((unsigned)time(NULL));

	vector<double> best(varCnt, 0);
	double minObjVal = (numeric_limits<double>::max)();

	//输入参数
	cout << "请输入目标函数的相关参数：" << endl;
	for (size_t i = 0; i < varCnt; i++)
	{
		double l, r, p;
		cout << "x" << i + 1 << "取值范围(以空格间开，e.g. -4 4)："; cin >> l >> r;
		cout << "x" << i + 1 << "精度(e.g. 0.001 or 1e-3)："; cin >> p;
		leftBound.push_back(l);
		rightBound.push_back(r);
		resolution.push_back(p);
	}
	cout << "请输入Nash Player的相关参数：" << endl;
	cout << "回合数(e.g. 60 or 6e1)："; cin >> RoundCnt;
	cout << "请输入GA的相关参数：" << endl;
	cout << "种群规模(e.g. 70 or 7e1)："; cin >> PopulationSize;
	cout << "迭代次数(e.g. 30 or 3e1)："; cin >> GenerationCnt;
	cout << "交叉概率(e.g. 0.65 or 6.5e-1)："; cin >> P_Cross;
	cout << "变异概率(e.g. 0.001 or 1e-3)："; cin >> P_Mutate;
	cout << "初始值[x0,y0](e.g. 3.14 -1.618 以空格间开)："; cin >> ans[0] >> ans[1];
	cout << "Calculating..." << endl << endl;

	//圆整
	for (size_t i = 0; i < varCnt; i++)
	{
		size_t chromoSegLen = ceil(log2l((rightBound[i] - leftBound[i]) / resolution[i]));
		resolution[i] = (rightBound[i] - leftBound[i]) / pow(2, chromoSegLen);
		varBitLen.push_back(chromoSegLen);
	}

	//参数信息写入文件
	fout << "当前测试参数如下：" << endl;
	fout << "目标函数：F(x,y)=((x-1)*y)^2 + (x*y-2)^4" << endl;
	for (size_t i = 0; i < varCnt; i++)
		fout << "x" << i + 1 << "取值区间：[" << leftBound[i] << "," << rightBound[i] << "]，基因串长度：" << varBitLen[i] << "，圆整后的实际精度：" << resolution[i] << endl;
	fout << "Nash 玩家数量：" << varCnt << endl;
	fout << "Nash 回合数：" << RoundCnt << endl;
	fout << "GA 种群规模：" << PopulationSize << endl;
	fout << "GA 迭代次数：" << GenerationCnt << endl;
	fout << "GA 交叉概率：" << P_Cross << endl;
	fout << "GA 变异概率：" << P_Mutate << endl << endl;

	cout << setw(6) << "Round" << setw(12) << "x" << setw(12) << "y" << setw(16) << "F(x,y)" << setw(8) << "Player" << endl;
	fout << setw(6) << "Round" << setw(12) << "x" << setw(12) << "y" << setw(16) << "F(x,y)" << setw(8) << "Player" << endl;

	//随机给定初始值
	//for (size_t i = 0; i < varCnt; i++)
	//	ans[i] = min((rand() % 100) / 100.0 * (rightBound[i] - leftBound[i]) + leftBound[i], rightBound[i]);

	//Nash GA 迭代
	for (size_t i = 0; i < RoundCnt; i++)
	{
		vector<vector<double>> local_var(varCnt, vector<double>(varCnt));
		
		double curObjVal = (*F)(ans);
		if (curObjVal < minObjVal)
		{
			minObjVal = curObjVal;
			best = ans;
		}

		cout << setw(6) << i << setw(12) << ans[0] << setw(12) << ans[1] << setw(16) << curObjVal << setw(8) << "prev" << endl;
		fout << setw(6) << i << setw(12) << ans[0] << setw(12) << ans[1] << setw(16) << curObjVal << setw(8) << "prev" << endl;

		for (size_t k = 0; k < varCnt; k++)
		{
			local_var[k] = ans;
			GA_Solver(F, local_var[k], k);

			curObjVal = (*F)(local_var[k]);
			if (curObjVal < minObjVal)
			{
				minObjVal = curObjVal;
				best = local_var[k];
			}

			cout << setw(6) << i << setw(12) << local_var[k][0] << setw(12) << local_var[k][1] << setw(16) << curObjVal << setw(8) << k << endl;
			fout << setw(6) << i << setw(12) << local_var[k][0] << setw(12) << local_var[k][1] << setw(16) << curObjVal << setw(8) << k << endl;
		}
		cout << endl;
		fout << endl;

		for (size_t k = 0; k < varCnt; k++)
			ans[k] = local_var[k][k];
	}

	//输出全局结果
	fout << endl << endl << "最优结果为：x=" << best[0] << "，y=" << best[1] << "，F(x,y)=" << minObjVal << endl;

	fout.close();
	return 0;
}
