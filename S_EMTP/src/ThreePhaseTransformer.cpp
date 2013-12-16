#include "ThreePhaseTransformer.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

ThreePhaseTransformer::ThreePhaseTransformer(int id,int firstNode,int lastNode,double Sbase,double Vprimary,double Vsecondary,double Frequency,double Xleakage,double r_lyw,double r_ayw)
{
	type=25;
	isUserDef = 0;
	need_NEC=1;
	nPort=7;
	for (int i=0;i<nPort;i++)
	{
		nodeNumberArray[i]= firstNode+i;
	}
	this->Sbase=Sbase;
	this->Vpbase=Vprimary;
	this->Vsbase=Vsecondary;
	this->Frequency=Frequency;
	this->Xleak=Xleakage;
	this->Rlyw=r_lyw;
	this->Alyw=r_ayw;
	this->Imag=0.01;
	this->Wbase=2*PI*Frequency;
	this->Ibase = Sbase/sqrt(3.0)/Vpbase;
	this->Rbase = Vpbase*Vpbase/Sbase;

	for (int i=0;i<nPort;i++)
	{
		nodeVoltage[i]=0;
		nodeVoltage_bak[i]=0;
		branchCurrent[i]=0;
		branchCurrent_bak[i]=0;
		nortonEquivalentCurrent[i]=0;
		nortonEquivalentCurrent_bak[i]=0;
		nortonEquivalentCurrent_1[i]=0;
		nortonEquivalentCurrent_2[i]=0;
	}
}

void ThreePhaseTransformer::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	readNodeVoltage(initialVoltageArray);
	for (int i=0;i<nPort;i++)
	{
		branchCurrent[i]= initialCurrentArray[ptr];
		ptr++;
	}

	for (int i=0;i<nPort;i++)
	{
		nortonEquivalentCurrent[i]=branchCurrent[i];
		for (int j=0;j<nPort;j++)
		{
			nortonEquivalentCurrent[i] += nortonEquivalentConductance(i+1,j+1)*nodeVoltage[j];
		}
		nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent[i];
	}
}

void ThreePhaseTransformer::readNodeVoltage(TVectorD& nodeVoltageArray)
{
	for(int i=0;i<nPort;i++)
	{
		if(nodeNumberArray[i]==0)
			nodeVoltage[i]=0;
		else
			nodeVoltage[i]=nodeVoltageArray[nodeNumberArray[i]];
	}
}

void ThreePhaseTransformer::calculateBranchCurrent()
{
	for (int i=0;i<nPort;i++)
	{
		branchCurrent_1[i]=branchCurrent[i];
		branchCurrent[i]=nortonEquivalentCurrent[i];
		for (int j=0;j<nPort;j++)
		{
			branchCurrent[i] += nortonEquivalentConductance(i+1,j+1)*nodeVoltage[j];
		}
	}
}

void ThreePhaseTransformer::formNodeNortonEquivalentCurrentArray(TVectorD& nodeNortonEquivalentCurrentArray)
{
	int N;
	for (int i=0;i<nPort;i++)
	{
		N = nodeNumberArray[i];
		nodeNortonEquivalentCurrentArray(N) = nodeNortonEquivalentCurrentArray(N)- nortonEquivalentCurrent[i];
	}
}

void ThreePhaseTransformer::formConductanceMatrix(TMatrixD& conductanceMatrix)
{
	int N,M;
	double tmpd;
	for (int i=0;i<nPort;i++)
	{
		N = nodeNumberArray[i];
		for (int j=0;j<nPort;j++)
		{
			M = nodeNumberArray[j];
			tmpd = conductanceMatrix(N,M);
			conductanceMatrix(N,M) = tmpd+this->nortonEquivalentConductance(i+1,j+1);
		}
	}
}


void ThreePhaseTransformer::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int i=0;i<nPort;i++)
	{
		branchCurrentMatrix_1[counter-1][ptr] = branchCurrent[i];
		ptr++;
	}
}

void ThreePhaseTransformer::interpolate(double ratio)
{
	for (int i=0;i<nPort;i++)
	{
		branchCurrent[i] = (1-ratio)*branchCurrent_1[i] + ratio*branchCurrent[i];
		nodeVoltage[i] = (1-ratio)*nodeVoltage_1[i] + ratio*nodeVoltage[i];
		nortonEquivalentCurrent[i] = (1-ratio)*nortonEquivalentCurrent_1[i] + ratio*nortonEquivalentCurrent[i];
	}
}

//更新开关处理过程中用于存储结果的变量
void ThreePhaseTransformer::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int i=0;i<nPort;i++)
		{
			branchCurrent_2[i] = branchCurrent_1[i];
			nodeVoltage_2[i] = nodeVoltage_1[i];
			nortonEquivalentCurrent_2[i]=nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<nPort;i++)
		{
			branchCurrent_1[i] = branchCurrent_2[i];
			nodeVoltage_1[i] = nodeVoltage_2[i];
			nortonEquivalentCurrent_1[i]=nortonEquivalentCurrent_2[i];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

// 存储内部变量，以便在迭代校正时恢复
void ThreePhaseTransformer::storeInternalVariables()
{
	for (int i=0;i<nPort;i++)
	{
		branchCurrent_bak[i] = branchCurrent[i];
		nodeVoltage_bak[i] = nodeVoltage[i];
		nortonEquivalentCurrent_bak[i]=nortonEquivalentCurrent[i];
	}
}
// 迭代校正时恢复内部变量
void ThreePhaseTransformer::restoreInternalVariables()
{
	for (int i=0;i<nPort;i++)
	{
		branchCurrent[i] = branchCurrent_bak[i];
		nodeVoltage[i] = nodeVoltage_bak[i];
		nortonEquivalentCurrent[i]=nortonEquivalentCurrent_bak[i];
	}
}

void ThreePhaseTransformer::calculateNortonEquivalentResistance(double time)
{
	double N_1_3_5=Vpbase/sqrt(3.0);
	double N_2_4_6=Vsbase;
	double L_w=1,A_w=1;
	double L_y=Rlyw;
	double A_y=Alyw;

	double L_1_6=0.5*L_w;
	double A_1_6 = A_w; 
	double L_13_14 = 0.5*L_y; 
	double A_13_14 = A_y;

	double r_R_1_13 = (L_1_6/A_1_6)/(L_13_14/A_13_14);
	double R_oc_3 = (Wbase*N_1_3_5*N_1_3_5*Imag*Ibase) /( Vpbase/sqrt(3.0));
	double R_13_14 = R_oc_3 / (2*r_R_1_13+2);
	double R_1_6 = R_13_14 * r_R_1_13;

	double P_1_6 = 1/R_1_6;
	double P_13_14 = 1/R_13_14;
	double R_st = (Wbase*N_1_3_5*N_1_3_5)/(Xleak*Rbase);
	double R_sigma1 = 2*R_st;
	double P_sigma1 =  1/R_sigma1;
	double P_sigma2 = P_sigma1;

	TMatrixD A;
	A.ResizeTo(1,17,1,6);
	A(1,1) = 1; A(7,1) = -1; A(13,1) = -1; A(15,1) = -1;
	A(1,2) = -1; A(2,2) = 1; A(7,2) = 1; A(8,2) = -1;
	A(3,3) = 1; A(9,3) = -1; A(13,3) = 1; A(14,3) = 1; A(16,3) = -1;
	A(3,4) = -1; A(4,4) = 1; A(9,4) = 1; A(10,4) = -1;
	A(5,5) = 1; A(11,5) = -1; A(14,5) = -1; A(17,5) = -1;
	A(5,6) = -1; A(6,6) = 1; A(11,6) = 1; A(12,6) = -1;

	TMatrixD P;
	P.ResizeTo(1,17,1,17);
	for (int i=1;i<=6;i++)
	{
		P(i,i) = P_1_6;
	}
	for (int i=7;i<=12;i++)
	{
		P(i,i) = P_sigma1;
	}
	for (int i=13;i<=14;i++)
	{
		P(i,i) = 0.5*P_13_14;
	}
	for (int i=15;i<=17;i++)
	{
		P(i,i) = P_sigma2;
	}

	TMatrixD I17;
	I17.ResizeTo(1,17,1,17);
	for (int i=1;i<=17;i++)
	{
		I17(i,i)=1;
	}

	TMatrixD At;
	At.ResizeTo(1,6,1,17);
	At.Transpose(A);
	TMatrixD Temp= At*P*A;
	TMatrixD Temp_inv = Temp;
	Temp_inv.Invert();
	TMatrixD M = (I17 - P*A * Temp_inv*At) * P;
	TMatrixD  M_ss= M;
	M_ss.ResizeTo(1,6,1,6);
	TMatrixD N_ss;
	N_ss.ResizeTo(1,6,1,6);
	for (int i=0;i<3;i++)
	{
		N_ss(2*i+1,2*i+1)=N_1_3_5;
		N_ss(2*i+2,2*i+2)=N_2_4_6;
	}
	TMatrixD I6;
	I6.ResizeTo(1,6,1,6);
	for (int i=1;i<=6;i++)
	{
		I6(i,i)=1;
	}

	TMatrixD MN_ss = M_ss*N_ss;
	TMatrixD MN_ss_inv=MN_ss.Invert();
	TMatrixD N_ss_inv=N_ss.Invert();
	TMatrixD Y_ss = (I6*MN_ss_inv) * 0.5*deltaT * (I6*N_ss_inv);
	TMatrixD Y;
	Y.ResizeTo(1,12,1,12);
	for(int i=1;i<=6;i++)
	{
		for (int j=1;j<=6;j++)
		{
			Y(i,j)= Y_ss(i,j);
		}
	}
	for(int i=7;i<=12;i++)
	{
		for (int j=1;j<=6;j++)
		{
			Y(i,j)= -Y_ss(i-6,j);
		}
	}
	for(int j=7;j<=12;j++)
	{
		for (int i=1;i<=6;i++)
		{
			Y(i,j)= -Y_ss(i,j-6);
		}
	}
	for(int i=7;i<=12;i++)
	{
		for (int j=7;j<=12;j++)
		{
			Y(i,j)= Y_ss(i-6,j-6);
		}
	}

	for (int i=1;i<=12;i++)
	{
		Y(7,i) = Y(7,i) + Y(9,i) + Y(11,i);//一次侧Y， 7&9&11连一起
		Y(4,i) = Y(4,i) + Y(8,i);// 二次侧Delta % 4-8连
		Y(6,i) = Y(6,i) + Y(10,i);//6-10连
		Y(2,i) = Y(2,i) + Y(12,i);//2-12连
	}
	for (int i=1;i<=12;i++)
	{
		Y(i,7) = Y(i,7) + Y(i,9) + Y(i,11);
		Y(i,4) = Y(i,4) + Y(i,8);
		Y(i,6) = Y(i,6) + Y(i,10);
		Y(i,2) = Y(i,2) + Y(i,12);
	}

	nortonEquivalentConductance.ResizeTo(1,12,1,12);
	nortonEquivalentConductance = Y;
	nortonEquivalentConductance.ResizeTo(1,nPort,1,nPort);
}

void ThreePhaseTransformer::calculateNortonEquivalentCurrent(double time)
{
	for (int i=0;i<nPort;i++)
	{
		nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent[i];
		nortonEquivalentCurrent[i]=branchCurrent[i];
		for (int j=0;j<nPort;j++)
		{
			nortonEquivalentCurrent[i] += nortonEquivalentConductance(i+1,j+1)*nodeVoltage[j];
		}
	}
}