#include "WoundInducGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                构造函数与析构函数                                    */
/************************************************************************/
WoundInducGenerator::WoundInducGenerator(int firstNode,int lastNode, int control,double vw)
{
	isUserDef = 0;
	need_NEC=1;
	//构造函数，a、b、c为节点编号
	nPort=6;
	nodeNumber=new int[nPort];
	type=18;

	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//用于初始化的参数
	w = 2*PI*50;//电角速度，有名值,f=60Hz
	w_pu=1;
	
	//dq0坐标系下的电阻参数,标幺值
	Rs = 0.048;
	Rr = 0.018;
	//dq0坐标系下的电抗参数,标幺值
	Xm = 3.8;
	Xls = 0.075;
	Xlr = 0.12;
	
	//定子侧电压、电流和阻抗的基值
	Vbase=690;
	Sbase=2e6;
	Vbase = Vbase/sqrt(3.0);//V
	Ibase = Sbase/3/Vbase;//A
	Zbase = Vbase/Ibase;
	
	//定转子绕组匝数比
	turnratio=0.3;

	this->control = control;

	//机械侧参数
	H = 3.0;
	//H = 0.25;
	D = 0;
	TL =-1.0;//输入转矩

	//风机参数
	p_air=1.225;
	rArea=3177;
	GR=90.5;
	beta=0;
	GE=0.979;
	vWind=vw;

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/
	//计算coff和T
	coff = 99.0/101.0;
	//coff=1;
	//coff = 0.8;
	T = 1.0/deltaT/w*(1+coff);

	//派克变换矩阵及其逆矩阵
	Ps.ResizeTo(1,2,1,3);
	invPs.ResizeTo(1,3,1,2);
	Pr.ResizeTo(1,2,1,3);
	invPr.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,4,1,4);
	BB.ResizeTo(1,4,1,4);
	////
	AA_ave.ResizeTo(1,4,1,4);
	AA_res.ResizeTo(1,4,1,4);
	inv_Aa.ResizeTo(1,4,1,4);
	////
	Gn1.ResizeTo(1,4,1,4);
	Gn2.ResizeTo(1,4,1,4);
	Gn3.ResizeTo(1,4,1,4);

	/************************************************************************/
	/*              各时步需要计算的物理量                                  */
	/************************************************************************/
	//以下各个物理量均为标幺值，先令其等于0，初始化时将赋初值
	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vsabc[i-1] = 0;
		Vsabc_1[i-1] = 0;
		Vsabc_2[i-1] = 0;
		Isabc[i-1] = 0;
		Isabc_1[i-1] = 0;
		Isabc_2[i-1] = 0;
	}
	//dq0坐标系下定子线圈的电压电流
	for (int i=1;i<=2;i++)
	{
		Vsdq0[i-1] = 0;
		Isdq0[i-1] = 0;
		Isdq0_his[i-1] = 0;
		Isdq0_his2[i-1] = 0;
		Isdq0_forcast[i-1] = 0;
	}

	//abc坐标系下转子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vrabc[i-1] = 0;
		Vrabc_1[i-1] = 0;
		Vrabc_2[i-1] = 0;
		Irabc[i-1] = 0;
		Irabc_1[i-1] = 0;
		Irabc_2[i-1] = 0;
	}
	//dq0坐标系下转子线圈的电压电流
	for (int i=1;i<=2;i++)
	{
		Vrdq0[i-1] = 0;
		Irdq0[i-1] = 0;
		Irdq0_his[i-1] = 0;
		Irdq0_his2[i-1] = 0;
		Irdq0_forcast[i-1] = 0;
	}

	/*  机械侧物理量  */
	if (control==1)
	{
		wr = 1.16715;//转子转速
		wr_his = 1.16715;//wr的历史值
	} 
	else
	{
		wr = 0;//转子转速
		wr_his = 0;//wr的历史值
	}
	curAngle_s= 0;//转子角,先赋为0，初始化时再重新赋值
	curAngle_r= 0;//转子角,先赋为0，初始化时再重新赋值


	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	for (int i=0;i<2;i++)
	{
		nortonEquivalentConductance[i] = 0;
		nortonEquivalentConductance_his[i] = 0;//等效导纳的历史值，用于更新节点导纳矩阵
	}
	
	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc坐标下的诺顿等效电流，有名值
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0坐标下的诺顿等效电流，标幺值
		nortonEquivalentCurrent_dq0_his[i-1] = 0;
	}

	/************************************************************************/
	/*              与EMTP接口所需的变量                                    */
	/************************************************************************/
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}

	/************************************************************************/
	/*              计算系数矩阵和诺等等值导纳                              */
	/************************************************************************/
	//calculateCoefficientMatrix();
}
WoundInducGenerator::~WoundInducGenerator()
{//析构函数

}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/

void WoundInducGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr)
{//初始化支路电压电流
	//初始化wr和curAngle
    //初始化**************************//
	//?????
	//**********************************
	
	T = 1.0/deltaT/w*(1+coff);
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
	}
	for (int i=0;i<3;i++)
	{
		Vrabc_1[i] =Vrabc[i];
	}

	for (int i=0;i<3;i++)
	{
		Isabc[i] = initialCurrentArray[ptr+i];
		Isabc_1[i] = Isabc[i];
	}
	ptr+=3;

	for (int i=0;i<3;i++)
	{
		Irabc[i] = initialCurrentArray[ptr+i];
		Irabc_1[i] = Irabc[i];
	}
	ptr+=3;

	calculateCoefficientMatrix();
	//解机械方程，初始化wr和curAngle

	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	for (int i=0;i<2;i++)
	{
		Vsdq0[i] = (Ps(i+1,1)*Vsabc[0]+Ps(i+1,2)*Vsabc[1]+Ps(i+1,3)*Vsabc[2])/Vbase;
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0_his2[i] = Vsdq0[i];
		Isdq0[i] = (Ps(i+1,1)*Isabc[0]+Ps(i+1,2)*Isabc[1]+Ps(i+1,3)*Isabc[2])/Ibase;
		Isdq0_his[i] = Isdq0[i];
		Isdq0_his2[i] = Isdq0[i];
		Isdq0_forcast[i] = Isdq0[i];
	}

	for (int i=0;i<2;i++)
	{
		Vrdq0[i] = (Pr(i+1,1)*Vrabc[0]+Pr(i+1,2)*Vrabc[1]+Pr(i+1,3)*Vrabc[2])/Vbase*turnratio;
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0_his2[i] = Vrdq0[i];
		Irdq0[i] = (Pr(i+1,1)*Irabc[0]+Pr(i+1,2)*Irabc[1]+Pr(i+1,3)*Irabc[2])/Ibase/turnratio;
		Irdq0_his[i] = Irdq0[i];
		Irdq0_his2[i] = Irdq0[i];
		Irdq0_forcast[i] = Irdq0[i];
	}

	//诺顿等值电流初始化
	//计算dq0坐标下的历史电流项
	for (int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			Gn1(i,1)*Isdq0_forcast[0]+Gn1(i,2)*Isdq0_forcast[1]
		+ Gn1(i,3)*Irdq0_forcast[0]+Gn1(i,4)*Irdq0_forcast[1]
		+Gn2(i,1)*Isdq0[0]+Gn2(i,2)*Isdq0[1]
		+Gn2(i,3)*Irdq0[0]+Gn2(i,4)*Irdq0[1]
		+ Gn3(i,1)*Vsdq0[0]+Gn3(i,2)*Vsdq0[1]
		+ Gn3(i,3)*Vrdq0[0]+Gn3(i,4)*Vrdq0[1];
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];
	}

	//计算park逆变换矩阵
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	//计算abc坐标下的历史电流项
	for (int i=1;i<=6;i++)
	{
		if (i<=3)
		{
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]
			+ invPs(i,2)*nortonEquivalentCurrent_dq0[1];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq0[2]
			+ invPr(i-3,2)*nortonEquivalentCurrent_dq0[3];
		    nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;
		}
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}
}

void WoundInducGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
{//从节点电压数组中读支路节点电压

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
	}

	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i]==0)
		{
			Vsabc[i]=0;
		}
		else
		{
			Vsabc[i]= nodeVoltageArray(nodeNumber[i]);
		}	
	}

	for (int i=0;i<3;i++)
	{
		Vrabc_1[i] =Vrabc[i];
	}


	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i+3]==0)
		{
			Vrabc[i]=0;
		}
		else
		{
			Vrabc[i]= nodeVoltageArray(nodeNumber[i+3]);
		}	
	}
}
void WoundInducGenerator::calculateBranchVoltage()
{//计算支路电压
}

#ifdef PREDICT
void WoundInducGenerator::calculateBranchCurrent()
{//计算支路电流
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i] = Vsabc[i]*nortonEquivalentConductance[0]+nortonEquivalentCurrent[i];
		Irabc_1[i] = Irabc[i];
		Irabc[i] = Vrabc[i]*nortonEquivalentConductance[1]+nortonEquivalentCurrent[i+3];
	}
}
#else
void WoundInducGenerator::calculateBranchCurrent()
{//计算支路电流
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i] = nortonEquivalentCurrent[i];
		for (int j=0;j<3;j++) {
			Isabc[i] += Vsabc[j]*Yne_abc(i+1,j+1);
			Isabc[i] += Vrabc[j]*Yne_abc(i+1,j+4);
		}

		Irabc_1[i] = Irabc[i];
		Irabc[i] = nortonEquivalentCurrent[i+3];
		for (int j=0;j<3;j++) {
			Irabc[i] += Vsabc[j]*Yne_abc(i+4,j+1);
			Irabc[i] += Vrabc[j]*Yne_abc(i+4,j+4);
		}
	}
}
#endif

void WoundInducGenerator::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	//if (control==2)
	//{
	//	//给定阶跃变化的风速
	//	if (time<1.4-deltaT/2)
	//	{
	//		vWind=14.8;
	//	}
	//	else if (time<1.6-deltaT/2)
	//	{
	//		vWind=13;
	//	}
	//	else if (time<1.8-deltaT/2)
	//	{
	//		vWind=11;
	//	}
	//	else if (time<2.2-deltaT/2)
	//	{
	//		vWind=16.8;
	//	}
	//	else
	//	{
	//		vWind=14.8;
	//	}
	//}

	calculateCoefficientMatrix();
	calculateDq0Results();//计算dq0坐标系下的结果，预报定子电流
	calculateOmega();

	//计算dq0坐标下的历史电流项
#ifdef PREDICT
	for (int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];
		nortonEquivalentCurrent_dq0[i-1] = 
			Gn1(i,1)*Isdq0_forcast[0]+Gn1(i,2)*Isdq0_forcast[1]
		+ Gn1(i,3)*Irdq0_forcast[0]+Gn1(i,4)*Irdq0_forcast[1]
		+ Gn2(i,1)*Isdq0[0]+Gn2(i,2)*Isdq0[1]
		+ Gn2(i,3)*Irdq0[0]+Gn2(i,4)*Irdq0[1]
		+ Gn3(i,1)*Vsdq0[0]+Gn3(i,2)*Vsdq0[1]
		+ Gn3(i,3)*Vrdq0[0]+Gn3(i,4)*Vrdq0[1];
	}
#else
	for (int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];
		nortonEquivalentCurrent_dq0[i-1] =
			  Gn2(i,1)*Isdq0[0] + Gn2(i,2)*Isdq0[1]
			+ Gn2(i,3)*Irdq0[0] + Gn2(i,4)*Irdq0[1]
			+ Gn3(i,1)*Vsdq0[0] + Gn3(i,2)*Vsdq0[1]
			+ Gn3(i,3)*Vrdq0[0] + Gn3(i,4)*Vrdq0[1];
	}
#endif

	//计算park逆变换矩阵
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//计算abc坐标下的历史电流项
	for (int i=1;i<=6;i++)
	{
		if (i<=3)
		{
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]
			+ invPs(i,2)*nortonEquivalentCurrent_dq0[1];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq0[2]
			+ invPr(i-3,2)*nortonEquivalentCurrent_dq0[3];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;
		}		
	}
}

void WoundInducGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//形成节点诺顿等效电流向量
	int N;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}

	// nodeNortonEquivalentCurrentArray.Print();
}
#ifdef PREDICT
void WoundInducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//形成节点导纳阵
	int N;
	double tmpd;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		if (i<3)
		{
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[0];
		} 
		else
		{
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[1];
		}
	}
}
#else
void WoundInducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//形成节点导纳阵
	int N1, N2;
	for (int i=0;i<6;i++)
	{
		N1 = nodeNumber[i];
		for (int j=0; j<6; j++)
		{
			N2 = nodeNumber[j];
			conductanceMatrix(N1,N2) += Yne_abc(i+1,j+1);
		}
	}
}
#endif

void WoundInducGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix(counter,ptr) = Isabc[0];
	branchCurrentMatrix(counter,ptr+1) = Isabc[1];
	branchCurrentMatrix(counter,ptr+2) = Isabc[2];
	ptr+=3;
	branchCurrentMatrix(counter,ptr) = Irabc[0];
	branchCurrentMatrix(counter,ptr+1) = Irabc[1];
	branchCurrentMatrix(counter,ptr+2) = Irabc[2];
	ptr+=3;
}

void WoundInducGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix_1[counter-1][ptr] = Isabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Isabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Isabc[2];
	ptr+=3;
	branchCurrentMatrix_1[counter-1][ptr] = Irabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Irabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Irabc[2];
	ptr+=3;
}
void WoundInducGenerator::interpolate(double ratio)
{//给定比值，对支路的电压电流进行插值

	//电压电流变量插值
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
	}

	for (int k=0;k<2;k++)
	{
		Vsdq0[k] = (1-ratio)*Vsdq0_his[k] + ratio*Vsdq0[k];
		Isdq0[k] = (1-ratio)*Isdq0_his[k] + ratio*Isdq0[k];
		Vrdq0[k] = (1-ratio)*Vrdq0_his[k] + ratio*Vrdq0[k];
		Irdq0[k] = (1-ratio)*Irdq0_his[k] + ratio*Irdq0[k];
	}

	for (int k=0;k<4;k++)
	{
		nortonEquivalentCurrent_dq0[k] = (1-ratio)*nortonEquivalentCurrent_dq0_his[k] + ratio*nortonEquivalentCurrent_dq0[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent[k] = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	}

	curAngle_r = (1-ratio)*curAngle_r_1 + ratio*curAngle_r;
	curAngle_s = (1-ratio)*curAngle_s_1 + ratio*curAngle_s;
	wr = (1-ratio)*wr_his + ratio*wr;

 	//诺顿等值电流插值
	//double temp1,temp2;
	//for (int k=0;k<6;k++)
	//{
	//	temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	//	temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

	//	nortonEquivalentCurrent_1[k] = temp1;
	//	nortonEquivalentCurrent[k] = temp2;
	//}
}

void WoundInducGenerator::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_2[i] = Vsabc_1[i];
			Isabc_2[i] = Isabc_1[i];
			Vrabc_2[i] = Vrabc_1[i];
			Irabc_2[i] = Irabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_2[i+3] = nortonEquivalentCurrent_1[i+3];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_1[i] = Vsabc_2[i];
			Isabc_1[i] = Isabc_2[i];
			Vrabc_1[i] = Vrabc_2[i];
			Irabc_1[i] = Irabc_2[i];
		}
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}
		break;
	case 3://诺顿等值电流存储：将_1变量的数值存入_2变量中
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

void WoundInducGenerator::updateConductanceMatrix(TMatrixD & conductanceMatrix)
{//转子转速变化时更新导纳矩阵
	int N;
	double tmpd;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		if (i<3)
		{
			conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his[0]+this->nortonEquivalentConductance[0];
		} 
		else
		{	
			conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his[1]+this->nortonEquivalentConductance[1];
		}
	}
}

/************************************************************************/
/*                发电机内部调用的函数                                  */
/************************************************************************/
//计算系数矩阵AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs4,Gr1...Gr3
//同时重算nortonEquivalentConductance
#ifdef PREDICT
void WoundInducGenerator::calculateCoefficientMatrix()
{
	//// 速度预测
	//double wr = 2*this->wr - this->wr_his;

	//计算AA矩阵，即伴随矩阵
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w_pu*(Xm+Xls);
	AA(1,3) = T*Xm;
	AA(1,4) = -w_pu*Xm;

	AA(2,1) = w_pu*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,3) = w_pu*Xm;
	AA(2,4) = T*Xm;

	AA(3,1) = T*Xm;
	AA(3,2) = (wr-w_pu)*Xm;
	AA(3,3) = Rr+T*(Xm+Xlr);
	AA(3,4) = (wr-w_pu)*(Xm+Xlr);

	AA(4,1) = (w_pu-wr)*Xm;
	AA(4,2) = T*Xm;
	AA(4,3) = (w_pu-wr)*(Xm+Xlr);
	AA(4,4) = Rr+T*(Xm+Xlr);

	//计算BB矩阵
	BB(1,1) = coff*Rs-T*(Xm+Xls);
	BB(1,2) = -coff*w_pu*(Xm+Xls);
	BB(1,3) = -T*Xm;
	BB(1,4) = -coff*w_pu*Xm;

	BB(2,1) = coff*w_pu*(Xm+Xls);
	BB(2,2) = coff*Rs-T*(Xm+Xls);
	BB(2,3) = coff*w_pu*Xm;
	BB(2,4) = -T*Xm;

	BB(3,1) = -T*Xm;
	BB(3,2) = coff*(wr-w_pu)*Xm;
	BB(3,3) = coff*Rr-T*(Xm+Xlr);
	BB(3,4) = coff*(wr-w_pu)*(Xm+Xlr);

	BB(4,1) = coff*(w_pu-wr)*Xm;
	BB(4,2) = -T*Xm;
	BB(4,3) = coff*(w_pu-wr)*(Xm+Xlr);
	BB(4,4) = coff*Rr-T*(Xm+Xlr);

	double Rsn=Rs+T*(Xm+Xls);
	double Rrn=Rr+T*(Xm+Xlr);
	AA_ave(1,1) = Rsn;
	AA_ave(2,2) = Rsn;
	AA_ave(3,3) = Rrn;
	AA_ave(4,4) = Rrn;
	for (int i=1;i<=4;i++)
	{
		for(int j=1;j<=4;j++)
		{
			if (i!=j)
			{
				AA_ave(i,j)=0;
			}
		}
	}

	AA_res = AA - AA_ave;

	//计算诺顿等效电导
	nortonEquivalentConductance[0] = 1/Rsn;
	nortonEquivalentConductance[1] = 1/Rrn;
	inv_Aa = AA_ave;
	inv_Aa.Invert();
	//将等效从标幺值换算到有名值
	nortonEquivalentConductance[0] /= Zbase;
	nortonEquivalentConductance[1] = nortonEquivalentConductance[1]/Zbase*turnratio*turnratio;

	//计算Gn1...Gn3
	Gn1 = (-1.0)*inv_Aa*AA_res;
	Gn2 = (-1.0)*inv_Aa*BB;
	Gn3 = coff*inv_Aa;
}
#else
void WoundInducGenerator::calculateCoefficientMatrix()
{
	//计算AA矩阵，即伴随矩阵
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w_pu*(Xm+Xls);
	AA(1,3) = T*Xm;
	AA(1,4) = -w_pu*Xm;

	AA(2,1) = w_pu*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,3) = w_pu*Xm;
	AA(2,4) = T*Xm;

	AA(3,1) = T*Xm;
	AA(3,2) = (wr-w_pu)*Xm;
	AA(3,3) = Rr+T*(Xm+Xlr);
	AA(3,4) = (wr-w_pu)*(Xm+Xlr);

	AA(4,1) = (w_pu-wr)*Xm;
	AA(4,2) = T*Xm;
	AA(4,3) = (w_pu-wr)*(Xm+Xlr);
	AA(4,4) = Rr+T*(Xm+Xlr);

	//计算BB矩阵
	BB(1,1) = coff*Rs-T*(Xm+Xls);
	BB(1,2) = -coff*w_pu*(Xm+Xls);
	BB(1,3) = -T*Xm;
	BB(1,4) = -coff*w_pu*Xm;

	BB(2,1) = coff*w_pu*(Xm+Xls);
	BB(2,2) = coff*Rs-T*(Xm+Xls);
	BB(2,3) = coff*w_pu*Xm;
	BB(2,4) = -T*Xm;

	BB(3,1) = -T*Xm;
	BB(3,2) = coff*(wr-w_pu)*Xm;
	BB(3,3) = coff*Rr-T*(Xm+Xlr);
	BB(3,4) = coff*(wr-w_pu)*(Xm+Xlr);

	BB(4,1) = coff*(w_pu-wr)*Xm;
	BB(4,2) = -T*Xm;
	BB(4,3) = coff*(w_pu-wr)*(Xm+Xlr);
	BB(4,4) = coff*Rr-T*(Xm+Xlr);

	Yne_dq.ResizeTo(1,4,1,4);
	Yne_abc.ResizeTo(1,6,1,6);

	Yne_dq = AA;
	Yne_dq.Invert();

	//计算Gn1...Gn3
	Gn2 = (-1.0)*Yne_dq*BB;
	Gn3 = coff*Yne_dq;
 }
#endif

//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
//为下一时步计算诺顿等效电流做准备
#ifdef PREDICT
void WoundInducGenerator::calculateDq0Results()
{
	double Vsabc_pu[3];
	double Isabc_pu[3];
	double Vrabc_pu[3];
	double Irabc_pu[3];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		Vsabc_pu[i] = Vsabc[i]/Vbase;
		Isabc_pu[i] = Isabc[i]/Ibase;
		Vrabc_pu[i] = Vrabc[i]/Vbase*turnratio;
		Irabc_pu[i] = Irabc[i]/Ibase/turnratio;
	}

	//计算park变换矩阵
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	//计算电压
	for (int i=0;i<2;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq0_his2[i] = Vrdq0_his[i];
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//计算电流
	for (int i=0;i<2;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq0_his2[i] = Irdq0_his[i];
		Irdq0_his[i] = Irdq0[i];
		Irdq0[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}

	//预报下一时刻电流
	for (int i=0;i<2;i++)
	{
		 Isdq0_forcast[i] = 1.3333333*Isdq0[i]+0.3333333*Isdq0_his[i]-0.6666666*Isdq0_his2[i];
		 Irdq0_forcast[i] = 1.3333333*Irdq0[i]+0.3333333*Irdq0_his[i]-0.6666666*Irdq0_his2[i];
		 //Isdq0_forcast[i] = 2*Isdq0[i] - Isdq0_his[i];
		 //Irdq0_forcast[i] = 2*Irdq0[i] - Irdq0_his[i];
		 //Isdq0_forcast[i] = 1.25*Isdq0[i]+0.5*Isdq0_his[i]-0.75*Isdq0_his2[i];
		 //Irdq0_forcast[i] = 1.25*Irdq0[i]+0.5*Irdq0_his[i]-0.75*Irdq0_his2[i];
	}
}
#else
void WoundInducGenerator::calculateDq0Results()
{
	double Vsabc_pu[3];
	double Isabc_pu[3];
	double Vrabc_pu[3];
	double Irabc_pu[3];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		Vsabc_pu[i] = Vsabc[i]/Vbase;
		Isabc_pu[i] = Isabc[i]/Ibase;
		Vrabc_pu[i] = Vrabc[i]/Vbase*turnratio;
		Irabc_pu[i] = Irabc[i]/Ibase/turnratio;
	}

	//计算park变换矩阵
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	//计算电压
	for (int i=0;i<2;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq0_his2[i] = Vrdq0_his[i];
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//计算电流
	for (int i=0;i<2;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq0_his2[i] = Irdq0_his[i];
		Irdq0_his[i] = Irdq0[i];
		Irdq0[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}

	//计算abc坐标下的节点导纳矩阵
	TMatrixD pt(1,4,1,6);
	TMatrixD ipt(1,6,1,4);

	for (int i=1;i<=2;i++) {
		for (int j=1;j<=3;j++) {
			pt(i,j) = Ps(i,j);
			pt(i+2,j+3) = Pr(i,j);
		}
	}

	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int i=1;i<=3;i++) {
		for (int j=1;j<=2;j++) {
			ipt(i,j) = invPs(i,j);
			ipt(i+3,j+2) = invPr(i,j);
		}
	}

	Yne_abc = ipt*Yne_dq*pt;

	for (int i=1; i<=6; i++)
		for (int j=1; j<=6; j++)
			Yne_abc(i,j) /= Zbase;
}
#endif

//求解机械方程，计算转子转速
//为下一时步求解系数矩阵、计算诺顿等值导纳做准备
void WoundInducGenerator::calculateOmega()
{
	//采用风力机拖动方式
	if (control==2)
	{
		//计算风力机转矩
		calculateTwind();
		//风力机转矩为拖动转矩
		TL=-Twind;
	}

	wr_his = wr;
	curAngle_r_1= curAngle_r;
	curAngle_s_1=curAngle_s;

	if (control==0||control==2)
	{
		Te_his=Te;
		Te=0.5*Xm*(Irdq0[0]*Isdq0[1]-Irdq0[1]*Isdq0[0]);
		wr = ((Te+Te_his)/2-TL)*deltaT/2/H+wr_his;
		curAngle_r = curAngle_r + deltaT*w*(w_pu-(wr_his+wr)/2);
		curAngle_s = curAngle_s + deltaT*w_pu*w;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(w_pu-wr)*w;
		curAngle_s = curAngle_s + deltaT*w_pu*w;
	}
}

//计算派克变换与反变换矩阵
void WoundInducGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void WoundInducGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}


// 平均模型相关，xuyin，20121224
// 存储内部变量，以便在迭代校正时恢复
void WoundInducGenerator::storeInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc_bak[k] = Vsabc[k];
		Isabc_bak[k] = Isabc[k];
		Vrabc_bak[k] = Vrabc[k];
		Irabc_bak[k] = Irabc[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq0_bak[k] = Vsdq0[k];
		Isdq0_bak[k] = Isdq0[k];
		Isdq0_his_bak[k] = Isdq0_his[k];
		Isdq0_his2_bak[k] = Isdq0_his2[k];
		Vrdq0_bak[k] = Vrdq0[k];
		Irdq0_bak[k] = Irdq0[k];
		Irdq0_his_bak[k] = Irdq0_his[k];
		Irdq0_his2_bak[k] = Irdq0_his2[k];
	}

	wr_bak = wr;
	wr_his_bak = wr_his;
	Te_bak = Te;
	Te_his_bak=Te_his;
	curAngle_s_bak = curAngle_s;
	curAngle_r_bak = curAngle_r;
}

// 迭代校正时恢复内部变量
void WoundInducGenerator::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc[k] = Vsabc_bak[k];
		Isabc[k] = Isabc_bak[k];
		Vrabc[k] = Vrabc_bak[k];
		Irabc[k] = Irabc_bak[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq0[k] = Vsdq0_bak[k];
		Isdq0[k] = Isdq0_bak[k];
		Isdq0_his[k] = Isdq0_his_bak[k];
		Isdq0_his2[k] = Isdq0_his2_bak[k];
		Vrdq0[k] = Vrdq0_bak[k];
		Irdq0[k] = Irdq0_bak[k];
		Irdq0_his[k] = Irdq0_his_bak[k];
		Irdq0_his2[k] = Irdq0_his2_bak[k];
	}

	wr = wr_bak;
	wr_his = wr_his_bak;
	Te=Te_bak;
	Te_his=Te_his_bak;
	curAngle_s = curAngle_s_bak;
	curAngle_r = curAngle_r_bak;
}

// 电机初始化专用函数, xuyin, 20121226
void WoundInducGenerator::initializeGen(double** GenInitialMatrix, int& ptr)
{
	// 除了第一列为时间t外，每个电机对应GenInitialMatrix的15列
	// 1-3列：Vsabc，4-6列：Vrabc，7-9列：Isabc，10-12列：Irabc，13-15列：Tm，Theta，wr
	// GenInitialMatrix共三行，分别对应当前时步以及前两个时步的值（按照时间顺序）

	// 电气量数值*1000
	for (int i=0; i<3; i++)
		for (int j=0; j<12; j++)
			GenInitialMatrix[i][ptr+j] *= 1000;

	// 读取时间
	double t = GenInitialMatrix[2][0]; //当前时间的坐标为第3行第1列

	// 初始化当前时刻的abc相电压电流
	for (int k=0; k<3; k++) {
		Vsabc[k] = GenInitialMatrix[2][ptr+k];
		Vrabc[k] = GenInitialMatrix[2][ptr+3+k];
		Isabc[k] = GenInitialMatrix[2][ptr+6+k];
		Irabc[k] = GenInitialMatrix[2][ptr+9+k];
	}

	// 机械转矩、转速和转子角
	TL = GenInitialMatrix[2][ptr+12];
	wr = GenInitialMatrix[2][ptr+14];
	wr_his = GenInitialMatrix[1][ptr+14];
	curAngle_s = w*w_pu*t;
	//curAngle_s = 0; // curAngle_s的初值可以任意选取
	double Theta = GenInitialMatrix[2][ptr+13];
	curAngle_r = curAngle_s - Theta;

	// 初始化电压电流的dq分量
	double dqMatrix[3][8];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		for (int j=ptr;j<ptr+3;j++)
			GenInitialMatrix[i][j] /= Vbase;

		for (int j=ptr+3;j<ptr+6;j++)
			GenInitialMatrix[i][j] /= (Vbase/turnratio);

		for (int j=ptr+6;j<ptr+9;j++)
			GenInitialMatrix[i][j] /= Ibase;

		for (int j=ptr+9;j<ptr+12;j++)
			GenInitialMatrix[i][j] /= (Ibase*turnratio);
	}

	// 计算dq分量
	for (int i=0; i<3; i++)
	{
		// 计算park变换角度
		double theta_s = w*GenInitialMatrix[i][0];
		double theta_r = theta_s - GenInitialMatrix[i][ptr+13];

		// 计算park变换矩阵
		parkTransMatrix(Ps,theta_s);
		parkTransMatrix(Pr,theta_r);

		// park变换
		// 定子电压
		for (int j=0; j<2; j++)
		{
			dqMatrix[i][j] = Ps(j+1,1)*GenInitialMatrix[i][ptr] + Ps(j+1,2)*GenInitialMatrix[i][ptr+1]
				+ Ps(j+1,3)*GenInitialMatrix[i][ptr+2];
		}
		// 转子电压
		for (int j=2; j<4; j++)
		{
			dqMatrix[i][j] = Pr(j-1,1)*GenInitialMatrix[i][ptr+3] + Pr(j-1,2)*GenInitialMatrix[i][ptr+4]
				+ Pr(j-1,3)*GenInitialMatrix[i][ptr+5];
		}
		// 定子电流
		for (int j=4; j<6; j++)
		{
			dqMatrix[i][j] = Ps(j-3,1)*GenInitialMatrix[i][ptr+6] + Ps(j-3,2)*GenInitialMatrix[i][ptr+7]
				+ Ps(j-3,3)*GenInitialMatrix[i][ptr+8];
		}
		// 转子电压
		for (int j=6; j<8; j++)
		{
			dqMatrix[i][j] = Pr(j-5,1)*GenInitialMatrix[i][ptr+9] + Pr(j-5,2)*GenInitialMatrix[i][ptr+10]
				+ Pr(j-5,3)*GenInitialMatrix[i][ptr+11];
		}
	}
	
	// dq分量赋值
	for (int k=0; k<2; k++) {
		Vsdq0[k] = dqMatrix[2][k];
		Isdq0[k] = dqMatrix[2][k+4];
		Isdq0_his[k] = dqMatrix[1][k+4];
		Isdq0_his2[k] = dqMatrix[0][k+4];
		Vrdq0[k] = dqMatrix[2][k+2];
		Irdq0[k] = dqMatrix[2][k+6];
		Irdq0_his[k] = dqMatrix[1][k+6];
		Irdq0_his2[k] = dqMatrix[0][k+6];
	}

	// 节点导纳矩阵计算
	calculateCoefficientMatrix();
	Te_his=0.5*Xm*(Irdq0_his2[0]*Isdq0_his2[1]-Irdq0_his2[1]*Isdq0_his2[0]);
	Te=0.5*Xm*(Irdq0_his[0]*Isdq0_his[1]-Irdq0_his[1]*Isdq0_his[0]);

#ifndef PREDICT
	TMatrixD pt(1,4,1,6);
	TMatrixD ipt(1,6,1,4);

	for (int i=1;i<=2;i++) {
		for (int j=1;j<=3;j++) {
			pt(i,j) = Ps(i,j);
			pt(i+2,j+3) = Pr(i,j);
		}
	}

	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int i=1;i<=3;i++) {
		for (int j=1;j<=2;j++) {
			ipt(i,j) = invPs(i,j);
			ipt(i+3,j+2) = invPr(i,j);
		}
	}

	Yne_abc = ipt*Yne_dq*pt;

	for (int i=1; i<=6; i++)
		for (int j=1; j<=6; j++)
			Yne_abc(i,j) /= Zbase;
#endif

	ptr += 15;
}

// 在风力机输入模式下，计算风力机输出转矩
void WoundInducGenerator::calculateTwind()
{
	TSR=2.237*vWind/(wr*w/GR);
	Cp=0.5*(TSR-0.022*beta*beta-5.6)*exp(-0.17*TSR);
	Pwind=0.5*Cp*rArea*(vWind*vWind*vWind)*p_air*GE/Sbase;
	Twind=Pwind/wr;
}