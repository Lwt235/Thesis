#include <iostream>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <set>
#include <algorithm>
#include <map>
#define LL long long
#define mod 998244353
using namespace std;
inline long long read() {
    long long x = 0, f = 1;
    char ch = getchar();
    while (ch < '0' || ch > '9')
    {
        if (ch == '-') f = -1;
        ch = getchar();
    }
    while (ch >= '0' && ch <= '9')
    {
        x = (x << 3) + (x << 1) + (ch ^ 48);
        ch = getchar();
    }
    return x * f;
}

const double trustLevel = 0.9;  //信任阈值
const double distrustLevel = 0.1; //不信任阈值
class individual {
public:
    string ID;
    vector<double> data;
    double weight;
    individual() { weight = 1.0; };
    individual(string ID, vector<double>data, double weight) :ID(ID), data(data), weight(weight) {};
};
template<typename T>
T dif(T a, T b) {
    return (a - b) * (a - b);
}

//信任曲线平滑度 a_smoothness>=1 值越小越光滑
const double a_smoothness = 2;

//0<=difPercent<=1

//快速幂 a^b
double qpow(double a, int b) {
    double res = 1;
    while (b) {
        if (b & 1) res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

//拟合基值
double funR(double x) {
    return x * (1 - qpow(x - 1, 2 * a_smoothness)) + (1 - x) * qpow(x, 2 * a_smoothness);
}
//拟合缩放
double funF(double x) {
    double X = 2 * a_smoothness * x;
    return X > 1 ? 0.5 * funR(1) : (X < 0 ? 0.5 * funR(0) : 0.5 * funR(X));
}
//拟合反转低值
double funG(double x) {
    return funF(-x - 1.0 / 6 * (1 - 1.0 / (2 * a_smoothness)));
}
//拟合反转高值
double funH(double x) {
    return funG(x + 2.0 / 3 * (1 - 1.0 / (2 * a_smoothness))) + 0.5;
}
//拟合双值合并
double funM(double x) {
    return x < -0.5 ? funH(x) : funG(x);
}
//信任度曲线拟合函数
double funAns(double x) {
    return funM(x - 1);
}

//与真值的差值百分比
double calcDif(double A, double X) {
    return 1 - 2.0 / (exp((A - X) / (sqrt(sqrt(qpow(X, 3))) + 1)) + exp((-A + X) / (sqrt(sqrt(qpow(X, 3))) + 1)));
}

double getCredibility(double difPercent) {
    //双曲正弦拟合
    /*double b = -1.0 / (2.0 * log(1.0 / a_credibility + sqrt(1.0 / a_credibility / a_credibility + 1.0)));
    return 0.5 + b * log((2 * difPercent - 1) / a_credibility + sqrt(((2 * difPercent - 1) / a_credibility) * ((2 * difPercent - 1) / a_credibility) + 1));*/

    //光滑曲线拟合
    return funAns(difPercent);
}

class Info {
public:
    double credibility;
    int trustCnt, distrustCnt;
};
map<string, Info> InfoSet;
template<typename T>
vector<T> evaluateX(vector<individual>& workerSet, const T& e) {
    if (workerSet.empty()) return {};
    vector<T>Res(workerSet[0].data.size());
    for (auto& worker : workerSet) if (!InfoSet.count(worker.ID)) InfoSet[worker.ID].credibility = 0.5;
    for (int t = 0; t < workerSet[0].data.size(); t++) {
        vector<individual>trustWorker;
        for (auto& worker : workerSet) {
            if (InfoSet[worker.ID].credibility >= trustLevel)
                trustWorker.emplace_back(worker);
        }
        if (trustWorker.empty()) trustWorker = workerSet;
        T difSum = 0;
        T X = 0;
        for (auto& worker : trustWorker) X += worker.data[t];
        X /= trustWorker.size();
        if (trustWorker.size() != 1) {
            for (auto& worker : trustWorker) difSum += dif(worker.data[t], X);
            T frontX = X;
            do {
                T weightSum = 0, average = 0;
                for (auto& worker : trustWorker) {
                    if (dif(worker.data[t], X) == 0) worker.weight = 1;
                    else worker.weight = log(difSum) - log(dif(worker.data[t], X));
                    weightSum += worker.weight;
                    average += worker.weight * worker.data[t];
                }
                frontX = X;
                X = average / weightSum;
            } while (dif(frontX, X) > e);
        }
        for (auto& worker : workerSet) {
            double dif = getCredibility(calcDif(worker.data[t], X)) - InfoSet[worker.ID].credibility;
            double trustPercent = 0, distrustPercent = 0;
            if (InfoSet[worker.ID].trustCnt + InfoSet[worker.ID].distrustCnt == 0)
                trustPercent = distrustPercent = 1;
            else {
                trustPercent = InfoSet[worker.ID].trustCnt * 1.0 / (InfoSet[worker.ID].trustCnt + InfoSet[worker.ID].distrustCnt);
                distrustPercent = 1 - trustPercent;
            }
            if (InfoSet[worker.ID].credibility + dif >= trustLevel) InfoSet[worker.ID].trustCnt++;
            else if (InfoSet[worker.ID].credibility + dif <= distrustLevel) InfoSet[worker.ID].distrustCnt++;
            worker.weight = InfoSet[worker.ID].credibility += dif * (dif > 0 ? trustPercent : distrustPercent);
        }
        cout << "第 " << t + 1 << " 组数据计算结束后信任度：" << endl;
        for (auto& item : InfoSet)
            cout << item.first << ": " << item.second.credibility << endl;
        for (auto& worker : workerSet) Res[t] = X;
    }
    return Res;
}
int main() {
    while (true) {
        cout << "请输入提供数据人员人数 n： ";
        int n; cin >> n; vector<individual>Test(n);
        cout << "请输入每位人员提供数据数目 m： ";
        int m; cin >> m;
        cout << "请以空格作为间隔，按：1.人员名称；2~m.数据内容 顺序输入：" << endl;
        for (int i = 0; i < n; i++) {
            cin >> Test[i].ID;
            for (int j = 0; j < m; j++)
            {
                double tmp; cin >> tmp;
                Test[i].data.emplace_back(tmp);
            }
        }
        cout << "请输入信任人员数量： ";
        cin >> n;
        if (n) cout << "请输入信任人员名称： ";
        string str;
        for (int i = 1; i <= n; i++) {
            cin >> str;
            InfoSet[str].credibility = 1;
            InfoSet[str].trustCnt++;
        }
        vector<double>Res = evaluateX(Test, 0.1);
        cout << "真值：";
        for (int i = 0; i < Res.size(); i++)
            printf("%lf ", Res[i]);
        cout << endl;
    }
}
/*

XiaoMin -10.05 2.01 1.56 2 1.5 1
XiaoLiang 15 30 25 0 20 20
XiaoHong 0.5 1.25 0.75 5 1.2 0.5
XiaoQiang 16 29 26 100 61 21
XiaoMei -21.25 1.3 0.9 1 3 0

*/