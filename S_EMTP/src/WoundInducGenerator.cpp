#include "WoundInducGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                构造函数与析构函数                                    */
/************************************************************************/
WoundInducGenerator::WoundInducGenerator(int id, int firstNode,int lastNode, int control,double Vw)
{
	this->id = id;
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
	Wbase = 2*PI*50;//电角速度，有名值,f=60Hz
	w=1;
	
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
	H = 3;
	D = 0;
	TL =-1.0;//输入转矩

	//风机参数
	p_air=1.225;
	rArea=3177;
	GR=90.5;
	beta=0;
	GE=0.979;
	vWind=Vw-0.1*(id-1);

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/

	//派克变换矩阵及其逆矩阵
	Ps.ResizeTo(1,2,1,3);
	invPs.ResizeTo(1,3,1,2);
	Pr.ResizeTo(1,2,1,3);
	invPr.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,4,1,4);
	BB.ResizeTo(1,4,1,4);
	AA_inv.ResizeTo(1,4,1,4);

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
		Isrc[i-1]=0;
		Isrc_1[i-1]=0;
		Isrc_2[i-1]=0;
	}
	//dq0坐标系下定子线圈的电压电流
	for (int i=1;i<=2;i++)
	{
		Vsdq[i-1] = 0;
		Isdq[i-1] = 0;
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
		Irrc[i-1]=0;
		Irrc_1[i-1]=0;
		Irrc_2[i-1]=0;
	}
	//dq0坐标系下转子线圈的电压电流
	for (int i=1;i<=2;i++)
	{
		Vrdq[i-1] = 0;
		Irdq[i-1] = 0;
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
	Te=0;
	Te_his=0;


	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	for (int i=0;i<2;i++)
	{
		nortonEquivalentConductance[i] = 0;
		nortonEquivalentConductance_pu[i]=0;
		RCnortonEquivalentConductance[i]=0;
	}

	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc坐标下的诺顿等效电流，有名值
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
		RCnortonEquivalentCurrent[i-1]=0;
		RCnortonEquivalentCurrent_1[i-1]=0;
		RCnortonEquivalentCurrent_2[i-1]=0;
	}
	for(int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq[i-1] = 0;//dq0坐标下的诺顿等效电流，标幺值
		nortonEquivalentCurrent_dq_his[i-1] = 0;
	}

	/************************************************************************/
	/*              与EMTP接口所需的变量                                    */
	/************************************************************************/
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}
#ifdef WIND_VELOCITY_DATA_INPUT
	WindVelocityCounter=0;
	WindVelocityCounter_bak=0;
#endif
}
WoundInducGenerator::~WoundInducGenerator()
{//析构函数

}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/

void WoundInducGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time){ptr+=6;}

void WoundInducGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
{//从节点电压数组中读支路节点电压

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
		Vsabc_his2[i] = Vsabc_his[i];
		Vsabc_his[i] = Vsabc[i];
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
		Vrabc_his2[i] = Vrabc_his[i];
		Vrabc_his[i] = Vrabc[i];
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

void WoundInducGenerator::calculateBranchCurrent()
{//计算支路电流
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
#ifdef S_MACHINE_WITH_COMPENSATE_RESISTANCE
		Isabc[i] = Vsabc[i]*nortonEquivalentConductance[0]+nortonEquivalentCurrent[i];
#endif
#ifdef S_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
		Isabc[i] = nortonEquivalentCurrent[i];
#endif
		Irabc_1[i] = Irabc[i];
#ifdef R_MACHINE_WITH_COMPENSATE_RESISTANCE
		Irabc[i] = Vrabc[i]*nortonEquivalentConductance[1]+nortonEquivalentCurrent[i+3];
#endif
#ifdef R_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
		Irabc[i] = nortonEquivalentCurrent[i+3];
#endif
		
		Isrc_1[i] = Isrc[i];
		Isrc[i] = Vsabc[i]*RCnortonEquivalentConductance[0]+RCnortonEquivalentCurrent[i];
		Irrc_1[i] = Irrc[i];
		Irrc[i] = Vrabc[i]*RCnortonEquivalentConductance[1]+RCnortonEquivalentCurrent[i+3];
	}
}

void WoundInducGenerator::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	
	calculateDqResults();//计算dq0坐标系下的结果
	calculateOmega();
	calculateCoefficientMatrix();

	////计算ABC坐标系下的电压插值预测
	double Vsdq_pre[2],Vrdq_pre[2];
	for (int i=0;i<2;i++)
	{
		Vsdq_pre[i] = Ps(i+1,1)*(1.25*Vsabc[0]+0.5*Vsabc_his[0]-0.75*Vsabc_his2[0])+Ps(i+1,2)*(1.25*Vsabc[1]+0.5*Vsabc_his[1]-0.75*Vsabc_his2[1])
			+Ps(i+1,3)*(1.25*Vsabc[2]+0.5*Vsabc_his[2]-0.75*Vsabc_his2[2]);
		Vrdq_pre[i] = Pr(i+1,1)*(1.25*Vrabc[0]+0.5*Vrabc_his[0]-0.75*Vrabc_his2[0])+Pr(i+1,2)*(1.25*Vrabc[1]+0.5*Vrabc_his[1]-0.75*Vrabc_his2[1])
			+Pr(i+1,3)*(1.25*Vrabc[2]+0.5*Vrabc_his[2]-0.75*Vrabc_his2[2]);
		Vsdq_pre[i] = Vsdq_pre[i]/Vbase;
		Vrdq_pre[i] = Vrdq_pre[i]/Vbase*turnratio;
	}

	TMatrixD CC=AA_inv*BB;
	//计算dq0坐标下的历史电流项
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq_his[i-1] = nortonEquivalentCurrent_dq[i-1];
		
		//计算电机支路电流
#ifdef S_MACHINE_WITHOUT_INT
		//不采用电压插值预测
		nortonEquivalentCurrent_dq[i-1] =-(CC(i,1)*Isdq[0]+CC(i,2)*Isdq[1]+CC(i,3)*Irdq[0]+CC(i,4)*Irdq[1])
		+ 2*(AA_inv(i,1)*Vsdq[0]+AA_inv(i,2)*Vsdq[1]+AA_inv(i,3)*Vrdq[0]+AA_inv(i,4)*Vrdq[1]);
#endif

#ifdef S_ABC_MACHINE_WITH_INT
		//采用ABC轴电压插值预测
		nortonEquivalentCurrent_dq[i-1] =-(CC(i,1)*Isdq[0]+CC(i,2)*Isdq[1]+CC(i,3)*Irdq[0]+CC(i,4)*Irdq[1])
			+ AA_inv(i,1)*Vsdq_pre[0]+AA_inv(i,2)*Vsdq_pre[1]+AA_inv(i,3)*Vrdq_pre[0]+AA_inv(i,4)*Vrdq_pre[1]
			+AA_inv(i,1)*Vsdq[0]+AA_inv(i,2)*Vsdq[1]+AA_inv(i,3)*Vrdq[0]+AA_inv(i,4)*Vrdq[1];
#endif

#ifdef S_DQ_MACHINE_WITH_INT
		////采用DQ轴电压插值预测
		nortonEquivalentCurrent_dq[i-1] =-(CC(i,1)*Isdq[0]+CC(i,2)*Isdq[1]+CC(i,3)*Irdq[0]+CC(i,4)*Irdq[1])
			+ AA_inv(i,1)*(1.25*Vsdq[0]+0.5*Vsdq_his[0]-0.75*Vsdq_his2[0])+AA_inv(i,2)*(1.25*Vsdq[1]+0.5*Vsdq_his[1]-0.75*Vsdq_his2[1])
			+AA_inv(i,3)*(1.25*Vrdq[0]+0.5*Vrdq_his[0]-0.75*Vrdq_his2[0])+AA_inv(i,4)*(1.25*Vrdq[1]+0.5*Vrdq_his[1]-0.75*Vrdq_his2[1])
			+AA_inv(i,1)*Vsdq[0]+AA_inv(i,2)*Vsdq[1]+AA_inv(i,3)*Vrdq[0]+AA_inv(i,4)*Vrdq[1];
#endif

		//加入DQ轴补偿电流
#ifdef S_DQ_COMPENSATE_WITHOUT_INT
		nortonEquivalentCurrent_dq[i-1] =nortonEquivalentCurrent_dq[i-1]-Vsdq[i-1]*nortonEquivalentConductance_pu[0];
#endif
#ifdef S_DQ_COMPENSATE_WITH_INT
		nortonEquivalentCurrent_dq[i-1] =nortonEquivalentCurrent_dq[i-1]-(1.25*Vsdq[i-1]+0.5*Vsdq_his[i-1]-0.75*Vsdq_his2[i-1])*nortonEquivalentConductance_pu[0];
		//nortonEquivalentCurrent_dq[i-1] =nortonEquivalentCurrent_dq[i-1]-(4.0/3.0*Vsdq[i-1]+1.0/3.0*Vsdq_his[i-1]-2.0/3.0*Vsdq_his2[i-1])*nortonEquivalentConductance_pu[0];
#endif

		nortonEquivalentCurrent_dq_his[i+1] = nortonEquivalentCurrent_dq[i+1];
		
		//计算电机支路电流
#ifdef R_MACHINE_WITHOUT_INT	//不采用电压插值预测
		nortonEquivalentCurrent_dq[i+1] =-(CC(i+2,1)*Isdq[0]+CC(i+2,2)*Isdq[1]+CC(i+2,3)*Irdq[0]+CC(i+2,4)*Irdq[1])
			+ 2*(AA_inv(i+2,1)*Vsdq[0]+AA_inv(i+2,2)*Vsdq[1]+AA_inv(i+2,3)*Vrdq[0]+AA_inv(i+2,4)*Vrdq[1]);
#endif
#ifdef R_ABC_MACHINE_WITH_INT		//采用ABC轴电压插值预测
		nortonEquivalentCurrent_dq[i+1] =-(CC(i+2,1)*Isdq[0]+CC(i+2,2)*Isdq[1]+CC(i+2,3)*Irdq[0]+CC(i+2,4)*Irdq[1])
			+ AA_inv(i+2,1)*Vsdq_pre[0]+AA_inv(i+2,2)*Vsdq_pre[1]+AA_inv(i+2,3)*Vrdq_pre[0]+AA_inv(i+2,4)*Vrdq_pre[1]
			+AA_inv(i+2,1)*Vsdq[0]+AA_inv(i+2,2)*Vsdq[1]+AA_inv(i+2,3)*Vrdq[0]+AA_inv(i+2,4)*Vrdq[1];
#endif
#ifdef R_DQ_MACHINE_WITH_INT 	 //采用DQ轴电压插值
		nortonEquivalentCurrent_dq[i+1] =-(CC(i+2,1)*Isdq[0]+CC(i+2,2)*Isdq[1]+CC(i+2,3)*Irdq[0]+CC(i+2,4)*Irdq[1])
			+ AA_inv(i+2,1)*(1.25*Vsdq[0]+0.5*Vsdq_his[0]-0.75*Vsdq_his2[0])+AA_inv(i+2,2)*(1.25*Vsdq[1]+0.5*Vsdq_his[1]-0.75*Vsdq_his2[1])
			+AA_inv(i+2,3)*(1.25*Vrdq[0]+0.5*Vrdq_his[0]-0.75*Vrdq_his2[0])+AA_inv(i+2,4)*(1.25*Vrdq[1]+0.5*Vrdq_his[1]-0.75*Vrdq_his2[1])
			+AA_inv(i+2,1)*Vsdq[0]+AA_inv(i+2,2)*Vsdq[1]+AA_inv(i+2,3)*Vrdq[0]+AA_inv(i+2,4)*Vrdq[1];
#endif

		//加入DQ轴补偿电流
#ifdef R_DQ_COMPENSATE_WITHOUT_INT
		nortonEquivalentCurrent_dq[i+1] =nortonEquivalentCurrent_dq[i+1]-Vrdq[i-1]*nortonEquivalentConductance_pu[1];
#endif
#ifdef R_DQ_COMPENSATE_WITH_INT
		nortonEquivalentCurrent_dq[i+1] =nortonEquivalentCurrent_dq[i+1]-(1.25*Vrdq[i-1]+0.5*Vrdq_his[i-1]-0.75*Vrdq_his2[i-1])*nortonEquivalentConductance_pu[1];
		//nortonEquivalentCurrent_dq[i+1] =nortonEquivalentCurrent_dq[i+1]-(4.0/3.0*Vrdq[i-1]+1.0/3.0*Vrdq_his[i-1]-2.0/3.0*Vrdq_his2[i-1])*nortonEquivalentConductance_pu[1];
#endif
	}

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
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq[0]+ invPs(i,2)*nortonEquivalentCurrent_dq[1];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;

#ifdef S_ABC_COMPENSATE_WITHOUT_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-Vsabc[i-1]*nortonEquivalentConductance[0];//ABC轴补偿电流
#endif
#ifdef S_ABC_COMPENSATE_WITH_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(1.25*Vsabc[i-1]+0.5*Vsabc_his[i-1]-0.75*Vsabc_his2[i-1])*nortonEquivalentConductance[0];//ABC轴补偿电流，带插值
			//nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(4.0/3.0*Vsabc[i-1]+1.0/3.0*Vsabc_his[i-1]-2.0/3.0*Vsabc_his2[i-1])*nortonEquivalentConductance[0];//ABC轴补偿电流，带插值
#endif
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq[2]+ invPr(i-3,2)*nortonEquivalentCurrent_dq[3];			
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;

#ifdef R_ABC_COMPENSATE_WITHOUT_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-Vrabc[i-4]*nortonEquivalentConductance[1];//ABC轴补偿电流
#endif
#ifdef R_ABC_COMPENSATE_WITH_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(1.25*Vrabc[i-4]+0.5*Vrabc_his[i-4]-0.75*Vrabc_his2[i-4])*nortonEquivalentConductance[1];//ABC轴补偿电流，带插值
			//nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(4.0/3.0*Vrabc[i-4]+1.0/3.0*Vrabc_his[i-4]-2.0/3.0*Vrabc_his2[i-4])*nortonEquivalentConductance[1];//ABC轴补偿电流，带插值
#endif
		}		
	}

	//计算RC支路的诺顿等值电流
	double Pc[2];
	Pc[0]=deltaT/2/Crc[0];
	Pc[1]=deltaT/2/Crc[1];
	for (int i=0;i<3;i++)
	{
		RCnortonEquivalentCurrent[i] = (Rrc[0]-Pc[0])/(Rrc[0]+Pc[0])*Isrc[i] - 1/(Rrc[0]+Pc[0])*Vsabc[i];
		RCnortonEquivalentCurrent[i+3] = (Rrc[1]-Pc[1])/(Rrc[1]+Pc[1])*Irrc[i] - 1/(Rrc[1]+Pc[1])*Vrabc[i];
	}
}

void WoundInducGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//形成节点诺顿等效电流向量
	int N;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) = nodeNortonEquivalentCurrentArray(N)- nortonEquivalentCurrent[i] - RCnortonEquivalentCurrent[i];
	}
}
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
#ifdef S_MACHINE_WITH_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[0]+RCnortonEquivalentConductance[0];
#endif
#ifdef S_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+RCnortonEquivalentConductance[0];
#endif
		} 
		else
		{
#ifdef R_MACHINE_WITH_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[1]+RCnortonEquivalentConductance[1];
#endif
#ifdef R_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+RCnortonEquivalentConductance[1];
#endif
		}
	}
}

void WoundInducGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix(counter,ptr) = Isabc[0]+Isrc[0];
	branchCurrentMatrix(counter,ptr+1) = Isabc[1]+Isrc[1];
	branchCurrentMatrix(counter,ptr+2) = Isabc[2]+Isrc[2];
	ptr+=3;
	branchCurrentMatrix(counter,ptr) = Irabc[0]+Irrc[0];
	branchCurrentMatrix(counter,ptr+1) = Irabc[1]+Irrc[1];
	branchCurrentMatrix(counter,ptr+2) = Irabc[2]+Irrc[2];
	ptr+=3;
}

void WoundInducGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix_1[counter-1][ptr] = Isabc[0]+Isrc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Isabc[1]+Isrc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Isabc[2]+Isrc[2];
	ptr+=3;
	branchCurrentMatrix_1[counter-1][ptr] = Irabc[0]+Irrc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Irabc[1]+Irrc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Irabc[2]+Irrc[2];
	ptr+=3;
}
void WoundInducGenerator::interpolate(double ratio)
{//给定比值，对支路的电压电流进行插值

	//电压电流变量插值
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Isrc[k] = (1-ratio)*Isrc_1[k] + ratio*Isrc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
		Irrc[k] = (1-ratio)*Irrc_1[k] + ratio*Irrc[k];
	}

	for (int k=0;k<2;k++)
	{
		Vsdq[k] = (1-ratio)*Vsdq_his[k] + ratio*Vsdq[k];
		Isdq[k] = (1-ratio)*Isdq_his[k] + ratio*Isdq[k];
		Vrdq[k] = (1-ratio)*Vrdq_his[k] + ratio*Vrdq[k];
		Irdq[k] = (1-ratio)*Irdq_his[k] + ratio*Irdq[k];
	}

	for (int k=0;k<4;k++)
	{
		nortonEquivalentCurrent_dq[k] = (1-ratio)*nortonEquivalentCurrent_dq_his[k] + ratio*nortonEquivalentCurrent_dq[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent[k] = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		RCnortonEquivalentCurrent[k] = (1-ratio)*RCnortonEquivalentCurrent_1[k] + ratio*RCnortonEquivalentCurrent[k];
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
			Isrc_2[i] = Isrc_1[i];
			Vrabc_2[i] = Vrabc_1[i];
			Irabc_2[i] = Irabc_1[i];
			Irrc_2[i] = Irrc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_2[i+3] = nortonEquivalentCurrent_1[i+3];
			RCnortonEquivalentCurrent_2[i] = RCnortonEquivalentCurrent_1[i];
			RCnortonEquivalentCurrent_2[i+3] = RCnortonEquivalentCurrent_1[i+3];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_1[i] = Vsabc_2[i];
			Isabc_1[i] = Isabc_2[i];
			Isrc_1[i] = Isrc_2[i];
			Vrabc_1[i] = Vrabc_2[i];
			Irabc_1[i] = Irabc_2[i];
			Irrc_1[i] = Irrc_2[i];
		}
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
			RCnortonEquivalentCurrent[i] = RCnortonEquivalentCurrent_1[i];
			RCnortonEquivalentCurrent_1[i] = RCnortonEquivalentCurrent_2[i];
		}
		break;
	case 3://诺顿等值电流存储：将_1变量的数值存入_2变量中
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			RCnortonEquivalentCurrent[i] = RCnortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

/************************************************************************/
/*                发电机内部调用的函数                                  */
/************************************************************************/
//计算系数矩阵AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs4,Gr1...Gr3
//同时重算nortonEquivalentConductance
void WoundInducGenerator::calculateCoefficientMatrix()
{

	//计算AA矩阵，即伴随矩阵
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w*(Xm+Xls);
	AA(1,3) = T*Xm;
	AA(1,4) = -w*Xm;

	AA(2,1) = w*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,3) = w*Xm;
	AA(2,4) = T*Xm;

	AA(3,1) = T*Xm;
	AA(3,2) = (wr-w)*Xm;
	AA(3,3) = Rr+T*(Xm+Xlr);
	AA(3,4) = (wr-w)*(Xm+Xlr);

	AA(4,1) = (w-wr)*Xm;
	AA(4,2) = T*Xm;
	AA(4,3) = (w-wr)*(Xm+Xlr);
	AA(4,4) = Rr+T*(Xm+Xlr);

	//计算BB矩阵
	BB(1,1) = Rs-T*(Xm+Xls);
	BB(1,2) = -w*(Xm+Xls);
	BB(1,3) = -T*Xm;
	BB(1,4) = -w*Xm;

	BB(2,1) = w*(Xm+Xls);
	BB(2,2) = Rs-T*(Xm+Xls);
	BB(2,3) = w*Xm;
	BB(2,4) = -T*Xm;

	BB(3,1) = -T*Xm;
	BB(3,2) = (wr-w)*Xm;
	BB(3,3) = Rr-T*(Xm+Xlr);
	BB(3,4) = (wr-w)*(Xm+Xlr);

	BB(4,1) = (w-wr)*Xm;
	BB(4,2) = -T*Xm;
	BB(4,3) = (w-wr)*(Xm+Xlr);
	BB(4,4) = Rr-T*(Xm+Xlr);

	AA_inv = AA;
	AA_inv.Invert();
}

//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
//为下一时步计算诺顿等效电流做准备
void WoundInducGenerator::calculateDqResults()
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
		Vsdq_his2[i] = Vsdq_his[i];
		Vsdq_his[i] = Vsdq[i];
		Vsdq[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq_his2[i] = Vrdq_his[i];
		Vrdq_his[i] = Vrdq[i];
		Vrdq[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//计算电流
	for (int i=0;i<2;i++)
	{
		Isdq_his2[i] = Isdq_his[i];
		Isdq_his[i] = Isdq[i];
		Isdq[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq_his2[i] = Irdq_his[i];
		Irdq_his[i] = Irdq[i];
		Irdq[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}
}
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
		Te=0.5*Xm*(Irdq[0]*Isdq[1]-Irdq[1]*Isdq[0]);
		wr = ((Te+Te_his)/2-TL)*deltaT/2/H+wr_his;
		curAngle_r = curAngle_r + deltaT*Wbase*(w-(wr_his+wr)/2);
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(w-wr)*Wbase;
		curAngle_s = curAngle_s + deltaT*w*Wbase;
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
		Vsabc_his_bak[k] = Vsabc_his[k];
		Vsabc_his2_bak[k] = Vsabc_his2[k];
		Isabc_bak[k] = Isabc[k];
		Isrc_bak[k] = Isrc[k];
		Vrabc_bak[k] = Vrabc[k];		
		Vrabc_his_bak[k] = Vrabc_his[k];
		Vrabc_his2_bak[k] = Vrabc_his2[k];
		Irabc_bak[k] = Irabc[k];
		Irrc_bak[k] = Irrc[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq_bak[k] = Vsdq[k];
		Vsdq_his_bak[k] = Vsdq_his[k];
		Vsdq_his2_bak[k] = Vsdq_his2[k];
		Isdq_bak[k] = Isdq[k];
		Isdq_his_bak[k] = Isdq_his[k];
		Isdq_his2_bak[k] = Isdq_his2[k];
		Vrdq_bak[k] = Vrdq[k];
		Vrdq_his_bak[k] = Vrdq_his[k];
		Vrdq_his2_bak[k] = Vrdq_his2[k];
		Irdq_bak[k] = Irdq[k];
		Irdq_his_bak[k] = Irdq_his[k];
		Irdq_his2_bak[k] = Irdq_his2[k];
	}

	wr_bak = wr;
	wr_his_bak = wr_his;
	Te_bak = Te;
	Te_his_bak=Te_his;
	curAngle_s_bak = curAngle_s;
	curAngle_r_bak = curAngle_r;
	WindVelocityCounter_bak = WindVelocityCounter;
}

// 迭代校正时恢复内部变量
void WoundInducGenerator::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc[k] = Vsabc_bak[k];
		Vsabc_his[k] = Vsabc_his_bak[k];
		Vsabc_his2[k] = Vsabc_his2_bak[k];
		Isabc[k] = Isabc_bak[k];
		Isrc[k] = Isrc_bak[k];
		Vrabc[k] = Vrabc_bak[k];
		Vrabc_his[k] = Vrabc_his_bak[k];
		Vrabc_his2[k] = Vrabc_his2_bak[k];
		Irabc[k] = Irabc_bak[k];
		Irrc[k] = Irrc_bak[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq[k] = Vsdq_bak[k];
		Vsdq_his[k] = Vsdq_his_bak[k];
		Vsdq_his2[k] = Vsdq_his2_bak[k];
		Isdq[k] = Isdq_bak[k];
		Isdq_his[k] = Isdq_his_bak[k];
		Isdq_his2[k] = Isdq_his2_bak[k];
		Vrdq[k] = Vrdq_bak[k];
		Vrdq_his[k] = Vrdq_his_bak[k];
		Vrdq_his2[k] = Vrdq_his2_bak[k];
		Irdq[k] = Irdq_bak[k];
		Irdq_his[k] = Irdq_his_bak[k];
		Irdq_his2[k] = Irdq_his2_bak[k];
	}

	wr = wr_bak;
	wr_his = wr_his_bak;
	Te = Te_bak;
	Te_his = Te_his_bak;
	curAngle_s = curAngle_s_bak;
	curAngle_r = curAngle_r_bak;
	WindVelocityCounter = WindVelocityCounter_bak;
}

void WoundInducGenerator::calculateNortonEquivalentResistance(double time)
{
	//计算特征电抗
	double Lds,Ldr,nortonEquivalentResistance[2];
	Lds=(Xls+Xm*Xlr/(Xm+Xlr))*Zbase/Wbase;
	Ldr=Xlr*Zbase/Wbase;
	//计算补偿导纳
	nortonEquivalentResistance[0]=2*Lds/deltaT;
	nortonEquivalentResistance[1]=2*Ldr/deltaT;
	nortonEquivalentConductance[0]=1/nortonEquivalentResistance[0];
	nortonEquivalentConductance[1]=1/nortonEquivalentResistance[1]*turnratio*turnratio;
	//将等效从有名值换算到标幺值
	nortonEquivalentConductance_pu[0] = nortonEquivalentConductance[0]*Zbase;
	nortonEquivalentConductance_pu[1] = nortonEquivalentConductance[1]*Zbase/turnratio/turnratio;
	
	T = 2.0/deltaT/Wbase;
	//计算RC支路参数
	Rrc[0]=20*Zbase;
	Rrc[1]=Rrc[0]/turnratio/turnratio;
	Crc[0]=10*deltaT/Rrc[0];
	Crc[1]=10*deltaT/Rrc[1];
	//Crc[0]=10*(5e-6)/Rrc[0];
	//Crc[1]=10*(5e-6)/Rrc[1];
	RCnortonEquivalentConductance[0]=1/(Rrc[0]+deltaT/2/Crc[0]);
	RCnortonEquivalentConductance[1]=1/(Rrc[1]+deltaT/2/Crc[1]);
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
		Vsabc_his[k] = GenInitialMatrix[1][ptr+k];
		Vsabc_his2[k] = GenInitialMatrix[0][ptr+k];
		Vrabc[k] = GenInitialMatrix[2][ptr+3+k];
		Vrabc_his[k] = GenInitialMatrix[1][ptr+3+k];
		Vrabc_his2[k] = GenInitialMatrix[0][ptr+3+k];
		Isabc[k] = GenInitialMatrix[2][ptr+6+k];
		Irabc[k] = GenInitialMatrix[2][ptr+9+k];
	}

	// 机械转矩、转速和转子角
	TL = GenInitialMatrix[2][ptr+12];
	wr = GenInitialMatrix[2][ptr+14];
	wr_his = GenInitialMatrix[1][ptr+14];
	 curAngle_s = w*Wbase*t;
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
		double theta_s = w*Wbase*GenInitialMatrix[i][0];
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
		Vsdq[k] = dqMatrix[2][k];
		Vsdq_his[k] = dqMatrix[1][k];
		Vsdq_his2[k] = dqMatrix[0][k];
		Isdq[k] = dqMatrix[2][k+4];
		Isdq_his[k] = dqMatrix[1][k+4];
		Isdq_his2[k] = dqMatrix[0][k+4];
		Vrdq[k] = dqMatrix[2][k+2];
		Vrdq_his[k] = dqMatrix[1][k+2];
		Vrdq_his2[k] = dqMatrix[0][k+2];
		Irdq[k] = dqMatrix[2][k+6];
		Irdq_his[k] = dqMatrix[1][k+6];
		Irdq_his2[k] = dqMatrix[0][k+6];
	}

	// 节点导纳矩阵计算
	calculateCoefficientMatrix();
	Te_his=0.5*Xm*(Irdq_his2[0]*Isdq_his2[1]-Irdq_his2[1]*Isdq_his2[0]);
	Te=0.5*Xm*(Irdq_his[0]*Isdq_his[1]-Irdq_his[1]*Isdq_his[0]);

	ptr += 15;
}

// 在风力机输入模式下，计算风力机输出转矩
void WoundInducGenerator::calculateTwind()
{
#ifdef WIND_VELOCITY_DATA_INPUT
	setWindVelocity();
#endif
	TSR=2.237*vWind/(wr*Wbase/GR);
	Cp=0.5*(TSR-0.022*beta*beta-5.6)*exp(-0.17*TSR);
	Pwind=0.5*Cp*rArea*(vWind*vWind*vWind)*p_air*GE/Sbase;
	Twind=Pwind/wr;
}

//保存电机转速, By Gao Haixiang
void WoundInducGenerator::saveMachineWr(double** machineWrMatrix, int& ptr, int counter)
{
	machineWrMatrix[counter-1][ptr] = wr;
	ptr++;
}

//	将每个时步的风速数据保存在每台风机中 ,By Gao Haixiang
void WoundInducGenerator::getWindVelocityData(double** VwMatrix, int rows, int& ptr)
{
	WindVelocityVector = new double[rows];
	for (int i=0;i<rows;i++)
	{
		WindVelocityVector[i] = VwMatrix[i][ptr];
	}
	ptr++;
}

// 设定风速, By Gao Haixiang
void WoundInducGenerator::setWindVelocity()
{
	vWind = WindVelocityVector[WindVelocityCounter];
	WindVelocityCounter++;
}