#include "PWMConverter.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                构造函数与析构函数                                    */
/************************************************************************/
// 构造函数
//电气节点顺序: a,b,c,dc1,dc2；控制节点顺序:a,b,c,tri
PWMConverter::PWMConverter(int firstACNode,int firstDCNode, double R, double L, double C, int firstCtrlNode, int subType)
{
	this->subType = subType;

	isUserDef = 0; // 是否参加插值？应用平均模型的时候不必要
	need_NEC=1; // 需要计算诺顿等值电流

	nPort = 5;
	nodeNumber = new int[nPort];
	type = 24; // 编号

	// 电气端口编号
	nodeNumber[0] = firstACNode;
	nodeNumber[1] = firstACNode + 1;
	nodeNumber[2] = firstACNode + 2;
	nodeNumber[3] = firstDCNode;
	nodeNumber[4] = firstDCNode + 1;

	for (int i=0; i<4; i++)
	{
		ctrlNodeNumber[i] = i+firstCtrlNode;
	}

	// 元件参数
	this->R = R + 0.001; //考虑开关导通电阻
	this->L = L;
	this->C = C;

	// 诺顿等值相关
	// 导纳矩阵
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Y11[i][j] = 0;

	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			Y12[i][j] = 0;

	for (int i=0; i<2; i++)
		for (int j=0; j<3; j++)
			Y21[i][j] = 0;

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			Y22[i][j] = 0;

	for (int i=0; i<5; i++)
		for (int j=0; j<5; j++)
			Yne[i][j] = 0;

	// 等值电流
	for (int k=0; k<3; k++)
		Ine1[k] = 0;

	for (int k=0; k<2; k++)
		Ine2[k] = 0;

	for (int k=0; k<5; k++)
		Ine[k] = 0;

	// 端口电压电流
	// 电压
	for (int k=0; k<3; k++)
		Ub1[k] = 0;

	for (int k=0; k<2; k++)
		Ub2[k] = 0;

	for (int k=0; k<5; k++)
		Ub[k] = 0;

	// 电流
	for (int k=0; k<3; k++)
		Ib1[k] = 0;

	for (int k=0; k<2; k++)
		Ib2[k] = 0;

	for (int k=0; k<5; k++)
		Ib[k] = 0;

	// 平均模型参数
	for (int k = 0; k < 3; k++) {
		tao_s1[k] = 0;
		tao_s2[k] = 0;
		tao_s[k] = 0;
	}

	for (int k = 0; k < 3; k++) {
		tao_s1_1[k] = 0;
		tao_s2_1[k] = 0;
	}

	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			D[i][j] = 0;
}

// 析构函数
PWMConverter::~PWMConverter()
{
}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/
//初始化支路电压电流
void PWMConverter::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	// 初始化支路电压和电流
	// 电压
	readNodeVoltage(initialVoltageArray);
	for (int k=0; k<3; k++) {
		Ub1[k] = Ub[k];
	}
	for (int k=0; k<2; k++) {
		Ub2[k] = Ub[k+3];
	}

	// 电流
	for (int k=0; k<5; k++) {
		Ib[k] = initialCurrentArray(ptr+k);
	}
	ptr += 5;
	for (int k=0; k<3; k++) {
		Ib1[k] = Ib[k];
	}
	for (int k=0; k<2; k++) {
		Ib2[k] = Ib[k+3];
	}
}

//从节点电压数组中读支路两节点电压
void PWMConverter::readNodeVoltage(TVectorD& nodeVoltageArray)
{
	for (int i=0;i<5;i++)
	{
		if (nodeNumber[i]==0)
		{
			Ub[i] = 0;
		}
		else
		{
			Ub[i] = nodeVoltageArray(nodeNumber[i]);
		}	
	}
}

// 计算支路电压
void PWMConverter::calculateBranchVoltage()
{
	for (int k=0; k<3; k++) {
		Ub1[k] = Ub[k];
	}
	for (int k=0; k<2; k++) {
		Ub2[k] = Ub[k+3];
	}
}

// 计算支路电流
void PWMConverter::calculateBranchCurrent()
{
	// Ib1 = Y11*Ub1 + Y12*Ub2 + Ine1
	for (int i=0; i<3; i++) {
		Ib1[i] = Ine1[i];
		for (int j=0; j<3; j++)
			Ib1[i] += Y11[i][j]*Ub1[j];
		for (int j=0; j<2; j++)
			Ib1[i] += Y12[i][j]*Ub2[j];
	}

	// Ib2 = Y21*Ub1 + Y22*Ub2 + Ine2
	for (int i=0; i<2; i++) {
		Ib2[i] = Ine2[i];
		for (int j=0; j<3; j++)
			Ib2[i] += Y21[i][j]*Ub1[j];
		for (int j=0; j<2; j++)
			Ib2[i] += Y22[i][j]*Ub2[j];
	}

	// Ib
	for (int k=0; k<3; k++) {
		Ib[k] = Ib1[k];
	}

	for (int k=0; k<2; k++) {
		Ib[k+3] = Ib2[k];
	}

	//// debug
	//for (int k=0; k<5; k++) {
	//	cout << Ib[k] << "  ";
	//}
	//cout << endl;
}

// 计算元件诺顿等值电流
void PWMConverter::calculateNortonEquivalentCurrent(double time)
{
	// Ine1 = Y11*Ub1 + Y12*Ub2 + (1-(deltaT/2)*(R/L))/(1+(deltaT/2)*(R/L))*Ib1
	double tmp = (1-(deltaT/2)*(R/L))/(1+(deltaT/2)*(R/L));
	for (int i=0; i<3; i++) {
		Ine1[i] = tmp*Ib1[i];
		for (int j=0; j<3; j++)
			Ine1[i] += Y11[i][j]*Ub1[j];
		for (int j=0; j<2; j++)
			Ine1[i] += Y12[i][j]*Ub2[j];
	}

	if ( subType == 1 )
	{
		// Ine2 = - DT*Ine1 - (2*C/deltaT)*Ub2 - DT*Ib1 - Ib2
		tmp = -2*C/deltaT;
		for (int i=0; i<2; i++) {
			Ine2[i] = tmp*Ub2[i] - Ib2[i];
			for (int j=0; j<3; j++)
				Ine2[i] -= D[j][i]*(Ine1[j]+Ib1[j]);
		}
	}
	else
	{
		// Ine2 = - DT*Ine1 - (2*C/deltaT)*M2*Ub2 - DT*Ib1 - Ib2
		// M2 = [1,-1;-1,1]
		tmp = -2*C/deltaT;
		Ine2[0] = tmp*(Ub2[0]-Ub2[1]);
		Ine2[1] = tmp*(Ub2[1]-Ub2[0]);
		for (int i=0; i<2; i++) {
			Ine2[i] -= Ib2[i];
			for (int j=0; j<3; j++)
				Ine2[i] -= D[j][i]*(Ine1[j]+Ib1[j]);
		}
	}

	// Ine
	for (int k=0; k<3; k++) {
		Ine[k] = Ine1[k];
	}

	for (int k=0; k<2; k++) {
		Ine[k+3] = Ine2[k];
	}
}

// 形成节点诺顿等效电流向量
void PWMConverter::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{
	int N;
	for (int k=0; k<5; k++)
	{
		N = nodeNumber[k];
		if ( N != 0 ) { 
			nodeNortonEquivalentCurrentArray(N) -= Ine[k];
		}
	}
}

// 形成节点导纳阵
void PWMConverter::formConductanceMatrix(TMatrixD &conductanceMatrix)
{
	// 没有考虑端口接地的情况
	int N1, N2;
	double tmpd;
	for (int i=0; i<5; i++)
	{
		N1 = nodeNumber[i];
		for (int j=0; j<5; j++)
		{
			N2 = nodeNumber[j];
			tmpd = conductanceMatrix(N1,N2);
			conductanceMatrix(N1,N2) = tmpd + Yne[i][j];
		}
	}
}

// 保存支路电流
void PWMConverter::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int k=0; k<5; k++) {
		branchCurrentMatrix_1[counter-1][ptr+k] = Ib[k];
	}
	ptr += 5;
}

// 计算节点导纳矩阵
void PWMConverter::calculateYne()
{
	// Y11
	Y11[0][0] = deltaT/(2*L)/(1+deltaT/2*R/L);
	Y11[1][1] = Y11[0][0];
	Y11[2][2] = Y11[0][0];

	// Y12
	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			Y12[i][j] = -Y11[0][0]*D[i][j];

	// Y21
	for (int i=0; i<2; i++)
		for (int j=0; j<3; j++)
			Y21[i][j] = Y12[j][i];

	// Y22
	Y22[0][0] = 2*C/deltaT;
	Y22[1][1] = Y22[0][0];
	if ( subType == 1 ) // xuyin,20121221,subType=1表示直流侧中点接地，subType=2则不接地
	{
		Y22[0][1] = 0;
		Y22[1][0] = 0;
	}
	else
	{
		Y22[0][1] = -Y22[0][0];
		Y22[1][0] = -Y22[0][0];
	}

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			for (int k=0; k<3; k++)
				Y22[i][j] += Y11[0][0]*D[k][i]*D[k][j];

	// 由Y11,Y12,Y21,Y22得到Yne
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Yne[i][j] = Y11[i][j];

	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			Yne[i][j+3] = Y12[i][j];

	for (int i=0; i<2; i++)
		for (int j=0; j<3; j++)
			Yne[i+3][j] = Y21[i][j];

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			Yne[i+3][j+3] = Y22[i][j];
}

#ifndef GENERAL
// 初始化开关时刻,version 1, 非通用版
void PWMConverter::initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr)
{
	// PWMInitialMatrix中每3列一组，对应一个PWM变流器，为a,b,c三相在一个开关周期内的开关脉冲信号
	// nStep表示每个开关周期包含的时步数，即为PWMInitialMatrix的行数减1

	// 初始化开关时刻tao_s1,tao_s2,tao_s
	for (int k=0; k<3; k++) {
		tao_s1[k] = 0.5;
		for (int i=0; i<=nStep/2; i++) {
			if (PWMInitialMatrix[i][ptr+k] == 1) {
				tao_s1[k] = i/(double)nStep;
				break;
			}
		}

		tao_s2[k] = 1.0;
		for (int i=nStep/2; i<=nStep; i++) {
			if (PWMInitialMatrix[i][ptr+k] == 0) {
				tao_s2[k] = i/(double)nStep;
				break;
			}
		}

		tao_s[k] = tao_s2[k] - tao_s1[k];
	}
	ptr += 3;

	for (int k=0; k<3; k++) {
		tao_s1_1[k] = tao_s1[k];
		tao_s2_1[k] = tao_s2[k];
	}

	// 计算矩阵D
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}

	//// debug
	//for (int k=0; k<3; k++) {
	//	cout<<tao_s1[k]<<"\t"<<tao_s2[k]<<"\t"<<tao_s[k]<<endl;
	//}

	// 计算诺顿等值导纳矩阵
	calculateYne();
}

#else
// 初始化开关时刻,version 2, 通用版, xuyin 20121211
void PWMConverter::initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr)
{
	// PWMInitialMatrix中每3列一组，对应一个PWM变流器，为a,b,c三相在一个开关周期内的开关脉冲信号
	// nStep表示每个开关周期包含的时步数，即为PWMInitialMatrix的行数减1

	// 初始化tao_s
	for (int k=0; k<3; k++) {
		tao_s[k] = 0.0;
		for (int i=1; i<=nStep; i++) {
			if ( PWMInitialMatrix[i-1][ptr+k] == 1 && PWMInitialMatrix[i][ptr+k] == 1 )
			{
					tao_s[k] += 1;
			}
			else if ( PWMInitialMatrix[i-1][ptr+k] == 0 && PWMInitialMatrix[i][ptr+k] == 0 )
			{
					// tao_s_new[k] += 0;
			}
			else
			{
				tao_s[k] += 0.5;
			}
		}
		tao_s[k] /= nStep;
	}
	ptr += 3;

	// 计算矩阵D
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}

	// 计算诺顿等值导纳矩阵
	calculateYne();
}
#endif

void PWMConverter::ctrlNodeNumber_PWM(int* ctrlNodeNumber, int k)
{
	for (int i=0; i<4; i++)
	{
		ctrlNodeNumber[4*k+i] = this->ctrlNodeNumber[i];
	}
}

// 预测开关时刻1：线性外插
void PWMConverter::predictSwithcingInstants()
{
	double tao_s1_new[3];
	double tao_s2_new[3];

	for (int k=0; k<3; k++) {
		tao_s1_new[k] = 2*tao_s1[k] - tao_s1_1[k];
		tao_s2_new[k] = 2*tao_s2[k] - tao_s2_1[k];
	}

	for (int k=0; k<3; k++) {
		tao_s1_1[k] = tao_s1[k];
		tao_s2_1[k] = tao_s2[k];
	}

	for (int k=0; k<3; k++) {
		tao_s1[k] = tao_s1_new[k];
		tao_s2[k] = tao_s2_new[k];
	}

	// 计算矩阵D
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}
	
	// 计算诺顿等值导纳矩阵
	calculateYne();
}

// 预测开关时刻2：迭代求解控制系统
int PWMConverter::predictSwithcingInstants(double ** PWMPredictMatrix, int nStep, int& ptr)
{
	// PWMPredictMatrix中每4列一组，对应一个PWM变流器
	// 前3列为a,b,c三相的调制波信号，第4列为载波信号
	// nStep表示每个开关周期包含的时步数，即为PWMPredictMatrix的行数减1

	// 预测tao_s1,tao_s2,tao_s
	double temp1, temp2;
	double tao_s1_pre[3];
	double tao_s2_pre[3];
	for (int k=0; k<3; k++) {
		tao_s1_pre[k] = 0.5;
		for (int i=0; i<=nStep/2; i++) {
			if (PWMPredictMatrix[i][ptr+k] >= PWMPredictMatrix[i][ptr+3]) {
				if (i > 0) {
					temp1 = fabs(PWMPredictMatrix[i][ptr+k]-PWMPredictMatrix[i][ptr+3]);
					temp2 = fabs(PWMPredictMatrix[i-1][ptr+k]-PWMPredictMatrix[i-1][ptr+3]);
					tao_s1_pre[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s1_pre[k] = i/(double)nStep;
				break;
			}
		}

		tao_s2_pre[k] = 1.0;
		for (int i=nStep/2; i<=nStep; i++) {
			if (PWMPredictMatrix[i][ptr+k] <= PWMPredictMatrix[i][ptr+3]) {
				if (i > nStep/2) {
					temp1 = fabs(PWMPredictMatrix[i][ptr+k]-PWMPredictMatrix[i][ptr+3]);
					temp2 = fabs(PWMPredictMatrix[i-1][ptr+k]-PWMPredictMatrix[i-1][ptr+3]);
					tao_s2_pre[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s2_pre[k] = i/(double)nStep;
				break;
			}
		}
	}
	ptr += 4;

	// 比较前后两次预测值
	double maxErr = 0.0;
	for (int k=0; k<3; k++) {
		if ( fabs(tao_s1[k]-tao_s1_pre[k]) > maxErr )
			maxErr = fabs(tao_s1[k]-tao_s1_pre[k]);
		if ( fabs(tao_s2[k]-tao_s2_pre[k]) > maxErr )
			maxErr = fabs(tao_s2[k]-tao_s2_pre[k]);
	}

	if ( maxErr < 0.005 )
	{
		return 0;
	}
	else
	{
		for (int k=0; k<3; k++) {
			tao_s1[k] = tao_s1_pre[k];
			tao_s2[k] = tao_s2_pre[k];
			tao_s[k] = tao_s2[k] - tao_s1[k];
		}

		// 计算矩阵D
		for (int k=0; k<3; k++) {
			D[k][0] = tao_s[k];
			D[k][1] = 1 - tao_s[k];
		}
	
		// 计算诺顿等值导纳矩阵
		calculateYne();

		return 1;
	}
}

#ifndef GENERAL
// 校正开关时刻version1,非通用版
// 开关周期起始时刻为载波最大值时刻
// 平均时间段长度需要恰好等于一个开关周期
int PWMConverter::correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr)
{
	// PWMCorrectMatrix中每4列一组，对应一个PWM变流器
	// 前3列为a,b,c三相的调制波信号，第4列为载波信号
	// nStep表示每个开关周期包含的时步数，即为PWMCorrectMatrix的行数减1

	// debug
	//for (int k=0; k<=nStep; k++){
	//	cout << k+1 << "\t";
	//	for (int j=0; j<4; j++)
	//		cout << PWMCorrectMatrix[k][j] << "\t";
	//	cout << endl;
	//}
	//cout << endl;

	// 计算求解得到的开关时刻,version1,非通用版
	// 开关周期起始时刻为载波最大值时刻
	// 平均时间段长度需要恰好等于一个开关周期
	double temp1, temp2;
	double tao_s1_new[3], tao_s2_new[3];
	for (int k=0; k<3; k++) {
		tao_s1_new[k] = 0.5;
		for (int i=0; i<=nStep/2; i++) {
			if (PWMCorrectMatrix[i][ptr+k] >= PWMCorrectMatrix[i][ptr+3]) {
				if (i > 0) {
					temp1 = fabs(PWMCorrectMatrix[i][ptr+k]-PWMCorrectMatrix[i][ptr+3]);
					temp2 = fabs(PWMCorrectMatrix[i-1][ptr+k]-PWMCorrectMatrix[i-1][ptr+3]);
					tao_s1_new[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s1_new[k] = i/(double)nStep;
				break;
			}
		}

		tao_s2_new[k] = 1.0;
		for (int i=nStep/2; i<=nStep; i++) {
			if (PWMCorrectMatrix[i][ptr+k] <= PWMCorrectMatrix[i][ptr+3]) {
				if (i > nStep/2) {
					temp1 = fabs(PWMCorrectMatrix[i][ptr+k]-PWMCorrectMatrix[i][ptr+3]);
					temp2 = fabs(PWMCorrectMatrix[i-1][ptr+k]-PWMCorrectMatrix[i-1][ptr+3]);
					tao_s2_new[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s2_new[k] = i/(double)nStep;
				break;
			}
		}
	}
	ptr += 4;

	//// debug
	//cout << "tao_s1:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s1[k] << "  ";
	//cout << endl;

	//cout << "tao_s2:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s2[k] << "  ";
	//cout << endl;

	//cout << "tao_s1_new:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s1_new[k] << "  ";
	//cout << endl;

	//cout << "tao_s2_new:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s2_new[k] << "  ";
	//cout << endl << endl;

	// 计算最大误差
	double maxError = 0;
	for (int k=0; k<3; k++) {
		if ( fabs(tao_s1[k]-tao_s1_new[k]) > maxError )
			maxError = fabs(tao_s1[k]-tao_s1_new[k]);
		if ( fabs(tao_s2[k]-tao_s2_new[k]) > maxError )
			maxError = fabs(tao_s2[k]-tao_s2_new[k]);
	}

	// 校正
	if ( maxError >= tol ) // 误差大于设定值时需要校正
	{
		for (int k=0; k<3; k++) {
			tao_s1[k] = tao_s1_new[k];
			tao_s2[k] = tao_s2_new[k];
			tao_s[k] = tao_s2[k] - tao_s1[k];
			D[k][0] = tao_s[k];
			D[k][1] = 1 - tao_s[k];
		}		

		calculateYne();

		return 1;
	}
	else
		return 0;
}

#else
// 校正开关时刻version2, 通用版, xuyin 20121211
int PWMConverter::correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr)
{
	// PWMCorrectMatrix中每4列一组，对应一个PWM变流器
	// 前3列为a,b,c三相的调制波信号，第4列为载波信号
	// nStep表示每个开关周期包含的时步数，即为PWMCorrectMatrix的行数减1

	//cout << "PWMCorrectMatrix:" << endl;
	//for (int i=0; i<=nStep; i++) {
	//	for (int j=0; j<4; j++)
	//		cout << PWMCorrectMatrix[i][j] << "   ";
	//	cout << endl;
	//}

	// 计算求解得到的开关时刻,version2,通用版
	double tao_s_new[3];
	for (int k=0; k<3; k++)
	{
		tao_s_new[k] = 0.0;
		for (int i=1; i<=nStep; i++)
		{
			if (PWMCorrectMatrix[i-1][ptr+k] >= PWMCorrectMatrix[i-1][ptr+3]
				&& PWMCorrectMatrix[i][ptr+k] >= PWMCorrectMatrix[i][ptr+3]) 
			{
				tao_s_new[k] += 1;
			}
			else if (PWMCorrectMatrix[i-1][ptr+k] < PWMCorrectMatrix[i-1][ptr+3]
				&& PWMCorrectMatrix[i][ptr+k] < PWMCorrectMatrix[i][ptr+3]) 
			{
				// tao_s_new[k] += 0;
			}
			else
			{
				double temp1 = PWMCorrectMatrix[i][ptr+k]-PWMCorrectMatrix[i][ptr+3];
				double temp2 = PWMCorrectMatrix[i-1][ptr+k]-PWMCorrectMatrix[i-1][ptr+3];
				tao_s_new[k] += 0.5*(1+(temp1+temp2)/(fabs(temp1)+fabs(temp2)));
			}
		}
		tao_s_new[k] /= nStep;
	}

	ptr += 4;

	// 计算最大误差
	double maxError = 0;
	for (int k=0; k<3; k++)
		if ( fabs(tao_s[k]-tao_s_new[k]) > maxError )
			maxError = fabs(tao_s[k]-tao_s_new[k]);

	if (ptr == 4)
		int mydebug = 0;

	// 校正
	if ( maxError >= tol ) // 误差大于设定值时需要校正
	{
		for (int k=0; k<3; k++) {
			tao_s[k] = tao_s[k] + 0.4*(tao_s_new[k]-tao_s[k]);
			D[k][0] = tao_s[k];
			D[k][1] = 1 - tao_s[k];
		}		

		calculateYne();

		return 1;
	}
	else
	{
		return 0;
	}
}
#endif

// 存储内部变量，以便在迭代校正时恢复
void PWMConverter::storeInternalVariables()
{
	// 对于PWM变流器，备份电压电流即可
	for (int k=0; k<3; k++) {
		Ub1_bak[k] = Ub1[k];
		Ib1_bak[k] = Ib1[k];
	}

	for (int k=0; k<2; k++) {
		Ub2_bak[k] = Ub2[k];
		Ib2_bak[k] = Ib2[k];
	}

	for (int k=0; k<5; k++) {
		Ub_bak[k] = Ub[k];
		Ib_bak[k] = Ib[k];
	}
}

// 迭代校正时恢复内部变量
void PWMConverter::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Ub1[k] = Ub1_bak[k];
		Ib1[k] = Ib1_bak[k];
	}

	for (int k=0; k<2; k++) {
		Ub2[k] = Ub2_bak[k];
		Ib2[k] = Ib2_bak[k];
	}

	for (int k=0; k<5; k++) {
		Ub[k] = Ub_bak[k];
		Ib[k] = Ib_bak[k];
	}
}

void  PWMConverter::correctYne(double* tao_s)
{
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}		

	calculateYne();
}