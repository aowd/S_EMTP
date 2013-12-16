#include "InducGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                构造函数与析构函数                                    */
/************************************************************************/
InducGenerator::InducGenerator(int a,int b,int c)
{
	isUserDef = 0;
	need_NEC=1;
	//构造函数，a、b、c为节点编号
	nPort=3;
	nodeNumber=new int[3];
	type=10;

	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//用于初始化的参数
	w = 2*PI*60;//电角速度，有名值,f=60Hz
	w_pu=1;
	
	//dq0坐标系下的电阻参数,标幺值
	Rs = 0.00621;
	Rr = 0.00025765;
	//dq0坐标系下的电抗参数,标幺值
	Xm = 0.34288;
	Xls = 0.15453;
	Xlr = 0.40075;
	
	//定子侧电压、电流和阻抗的基值
	Vbase=10e3;
	double Sbase=5e6;
	Vbase = (Vbase/sqrt(3.0))*sqrt(2.0);//V
	Ibase = Sbase/1.5/Vbase;//A
	Zbase = Vbase/Ibase;
	
	//定转子绕组匝数比
	turnratio=1;

	control=0;//0为转矩控制，1为转速控制

	//机械侧参数
	H = 1.0;
	D = 0;
	Tm =0.8;//输入转矩

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/
	//计算coff和T
	//coff = 99/101;
	coff =39.0/41.0;

	//派克变换矩阵及其逆矩阵
	Ps.ResizeTo(1,2,1,3);
	invPs.ResizeTo(1,3,1,2);
	Pr.ResizeTo(1,2,1,3);
	invPr.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,4,1,4);
	BB.ResizeTo(1,4,1,4);
	////
	Rss.ResizeTo(1,2,1,2);
	Rsr.ResizeTo(1,2,1,2);
	Rrs.ResizeTo(1,2,1,2);
	Rrr.ResizeTo(1,2,1,2);

	////
	Bss.ResizeTo(1,2,1,2);
	Bsr.ResizeTo(1,2,1,2);
	Brs.ResizeTo(1,2,1,2);
	Brr.ResizeTo(1,2,1,2);

	////处理系数矩阵，使得转换到abc坐标系时能解耦
	R_s.ResizeTo(1,2,1,2);
	R_ave.ResizeTo(1,2,1,2);
	R_res.ResizeTo(1,2,1,2);
	inv_Rs.ResizeTo(1,2,1,2);
	inv_Rr.ResizeTo(1,2,1,2);

	////
	Gs1.ResizeTo(1,2,1,2);
	Gs2.ResizeTo(1,2,1,2);
	Gs3.ResizeTo(1,2,1,2);
	Gs4.ResizeTo(1,2,1,2);

	////
	Gr1.ResizeTo(1,2,1,2);
	Gr2.ResizeTo(1,2,1,2);
	Gr3.ResizeTo(1,2,1,2);

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
	wr = 0;//转子转速
	wr_his = 0;//wr的历史值
	curAngle_s= 0;//转子角,先赋为0，初始化时再重新赋值
	curAngle_r= 0;//转子角,先赋为0，初始化时再重新赋值


	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	nortonEquivalentConductance = 0;
	nortonEquivalentConductance_his = 0;//等效导纳的历史值，用于更新节点导纳矩阵

	for(int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc坐标下的诺顿等效电流，有名值
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0坐标下的诺顿等效电流，标幺值
	}

	/************************************************************************/
	/*              与EMTP接口所需的变量                                    */
	/************************************************************************/
	nodeNumber[0] = a;
	nodeNumber[1] = b;
	nodeNumber[2] = c;

	/************************************************************************/
	/*              计算系数矩阵和诺等等值导纳                              */
	/************************************************************************/
	//calculateCoefficientMatrix();
}
InducGenerator::~InducGenerator()
{//析构函数

}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/

void InducGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time)
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
		Isabc[i] = initialCurrentArray[ptr+i];
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
		Vrdq0[i] = (Pr(i+1,1)*Vrabc[0]+Pr(i+1,2)*Vrabc[1]+Pr(i+1,3)*Vrabc[2])/Vbase;
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0_his2[i] = Vrdq0[i];
		Irdq0[i] = (Pr(i+1,1)*Irabc[0]+Pr(i+1,2)*Irabc[1]+Pr(i+1,3)*Irabc[2])/Ibase;
		Irdq0_his[i] = Irdq0[i];
		Irdq0_his2[i] = Irdq0[i];
		Irdq0_forcast[i] = Irdq0[i];
	}

	//诺顿等值电流初始化
	//计算dq0坐标下的历史电流项
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			Gs1(i,1)*Isdq0_forcast[0]+Gs1(i,2)*Isdq0_forcast[1]
		+ Gs2(i,1)*Isdq0[0]+Gs2(i,2)*Isdq0[1]
		+ Gs3(i,1)*Irdq0[0]+Gs3(i,2)*Irdq0[1]
		+ Gs4(i,1)*Vsdq0[0]+Gs4(i,2)*Vsdq0[1];
	}

	//计算park逆变换矩阵
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int k=0;k<3;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//计算abc坐标下的历史电流项
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]
		+ invPs(i,2)*nortonEquivalentCurrent_dq0[1];
		nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
	}

}

void InducGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
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
}
void InducGenerator::calculateBranchVoltage()
{//计算支路电压
}

void InducGenerator::calculateBranchCurrent()                                                          
{//计算支路电流
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i] = Vsabc[i]*nortonEquivalentConductance+nortonEquivalentCurrent[i];
	}
}

void InducGenerator::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	calculateCoefficientMatrix();
	calculateDq0Results();//计算dq0坐标系下的结果，预报定子电流
	calculateOmega();

	//计算dq0坐标下的历史电流项
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			Gs1(i,1)*Isdq0_forcast[0]+Gs1(i,2)*Isdq0_forcast[1]
		+ Gs2(i,1)*Isdq0[0]+Gs2(i,2)*Isdq0[1]
		+ Gs3(i,1)*Irdq0[0]+Gs3(i,2)*Irdq0[1]
		+ Gs4(i,1)*Vsdq0[0]+Gs4(i,2)*Vsdq0[1];
	}

	//计算park逆变换矩阵
	invParkTransMatrix(invPs,curAngle_s);

	for (int k=0;k<3;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//计算abc坐标下的历史电流项
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]
		+ invPs(i,2)*nortonEquivalentCurrent_dq0[1];
		nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
	}
}

void InducGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//形成节点诺顿等效电流向量
	int N;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}
}
void InducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//形成节点导纳阵
	int N;
	double tmpd;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance;
	}
}

void InducGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix(counter,ptr) = Isabc[0];
	branchCurrentMatrix(counter,ptr+1) = Isabc[1];
	branchCurrentMatrix(counter,ptr+2) = Isabc[2];
	ptr+=3;
}

void InducGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix_1[counter-1][ptr] = Isabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Isabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Isabc[2];
	ptr+=3;
}
void InducGenerator::interpolate(double ratio)
{//给定比值，对支路的电压电流进行插值

	//电压电流变量插值
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
	}

	//诺顿等值电流插值
	double temp1,temp2;
	for (int k=0;k<3;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void InducGenerator::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_2[i] = Vsabc_1[i];
			Isabc_2[i] = Isabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_1[i] = Vsabc_2[i];
			Isabc_1[i] = Isabc_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://诺顿等值电流存储：将_1变量的数值存入_2变量中
		for (int i=0;i<3;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

void InducGenerator::updateConductanceMatrix(TMatrixD & conductanceMatrix)
{//转子转速变化时更新导纳矩阵
	int N;
	double tmpd;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his+this->nortonEquivalentConductance;
	}
}

/************************************************************************/
/*                发电机内部调用的函数                                  */
/************************************************************************/
//计算系数矩阵AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs4,Gr1...Gr3
//同时重算nortonEquivalentConductance
void InducGenerator::calculateCoefficientMatrix()
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

	//计算矩阵Rss,Rsr,Rrs,Rrr,Bss,Bsr,Brs,Brr
	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rss(i,j) = AA(i,j);
			Bss(i,j) = BB(i,j);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rsr(i,j) = AA(i,j+2);
			Bsr(i,j) = BB(i,j+2);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rrs(i,j) = AA(i+2,j);
			Brs(i,j) = BB(i+2,j);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rrr(i,j) = AA(i+2,j+2);
			Brr(i,j) = BB(i+2,j+2);
		}
	}

	//计算Rs,R_ave,R_res,inv_Rs,inv_Rr
	inv_Rr = Rrr;
	inv_Rr.Invert();

	R_s = Rss-Rsr*inv_Rr*Rrs;

//	double resistance = (R_s(1,1)+R_s(2,2))/2.0;
	double resistance=1/0.047619047619048;
	R_ave(1,1) = resistance;
	R_ave(1,2) = 0;
	R_ave(2,1) = 0;
	R_ave(2,2) = resistance;

	R_res = R_s - R_ave;

	//计算诺顿等效电导
	nortonEquivalentConductance = 1/resistance;
	inv_Rs = R_ave;
	inv_Rs.Invert();
	nortonEquivalentConductance /= Zbase;

	//计算Gs1...Gs4
	Gs1 = (-1.0)*inv_Rs*R_res;
	Gs2 = (-1.0)*inv_Rs*Bss+inv_Rs*Rsr*inv_Rr*Brs;
	Gs3 = (-1.0)*inv_Rs*Bsr+inv_Rs*Rsr*inv_Rr*Brr;
	Gs4 = coff*inv_Rs;

	//计算Gr1...Gr3
	Gr1 = (-1.0)*inv_Rr*Rrs;
	Gr2 = (-1.0)*inv_Rr*Brs;
	Gr3 = (-1.0)*inv_Rr*Brr;
}

//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
//为下一时步计算诺顿等效电流做准备
void InducGenerator::calculateDq0Results()
{
	double Vsabc_pu[3];
	double Isabc_pu[3];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		Vsabc_pu[i] = Vsabc[i]/Vbase;
		Isabc_pu[i] = Isabc[i]/Ibase;
	}

	//计算park变换矩阵
	parkTransMatrix(Ps,curAngle_s);

	//计算定子电压
	for (int i=0;i<2;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
	}

	//计算定子电流
	for (int i=0;i<2;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
	}

	//预报下一时刻定子电流??????????????????????????????
	for (int i=0;i<2;i++)
	{
		Isdq0_forcast[i] = 1.3333333*Isdq0[i]+0.3333333*Isdq0_his[i]-0.6666666*Isdq0_his2[i];
		//Isdq0_forcast[i] = 1.25*Isdq0[i]+0.25*Isdq0_his[i]-0.5*Isdq0_his2[i];
		//Isdq0_forcast[i] = Isdq0[i];
	}

	//计算转子电压
	Vrdq0[0] = 0;
	Vrdq0[1] = 0;

	//计算转子电流
	for (int i=0;i<2;i++)
	{
		Irdq0_his[i] = Irdq0[i];
	}
	for (int i=1;i<=2;i++)
	{
		Irdq0[i-1] = Gr1(i,1)*Isdq0[0] + Gr1(i,2)*Isdq0[1]
		+ Gr2(i,1)*Isdq0_his[0] + Gr2(i,2)*Isdq0_his[1]
		+ Gr3(i,1)*Irdq0_his[0] + Gr3(i,2)*Irdq0_his[1];
	}
}

//求解机械方程，计算转子转速
//为下一时步求解系数矩阵、计算诺顿等值导纳做准备
void InducGenerator::calculateOmega()
{
	wr_his = wr;
	cout<<wr<<endl;
	if (control==0)
	{
		wr = (Tm-Xm*(Irdq0[0]*Isdq0[1]-Irdq0[1]*Isdq0[0]))*deltaT/2/H+wr_his;
		curAngle_r = curAngle_r + deltaT*w*(1-(wr_his+wr)/2);
		curAngle_s=curAngle_s+ deltaT*w;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(1-wr)*w;
		curAngle_s=curAngle_s+ deltaT*w;
	}

	if (curAngle_s>2*PI)
	{
		curAngle_s -= 2*PI;
	}
	if (curAngle_r>2*PI)
	{
		curAngle_r -= 2*PI;
	}
	if (curAngle_r<0)
	{
		curAngle_r += 2*PI;
	}
}

//计算派克变换与反变换矩阵
void InducGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void InducGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}
