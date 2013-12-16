#include "SimpleBranch.h"
#include <iostream>
using namespace std;

SimpleBranch::SimpleBranch()
{
	nodeNumber=new int[2];
	nPort=2;
	for(int i=0;i<2;i++)
	{
		nodeNumber[i] = 0;
		nodeVoltage[i] = 0;
	}
	branchVoltage = 0;
	branchVoltage_1 = 0;
	branchCurrent = 0;
	branchCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent_2 = 0;
	nortonEquivalentResistance = 0;
}

bool SimpleBranch::isSeries()
{
	int from = nodeNumber[0];
	int to = nodeNumber[1];
	int series = 1;

	if(to==0) {to=from;series=0;}
	if(from==0) {from=to;series=0;}
	if(to==0) cerr<<"Both Nodes Zero!!"<<endl;

	return series;
}

void SimpleBranch::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();
	branchCurrent = initialCurrentArray[ptr];
	ptr++;

	// edited by xuyin, 20121208
	branchCurrent_1 = branchCurrent;
	branchVoltage_1 = branchVoltage;
	nortonEquivalentCurrent = branchCurrent - branchVoltage/nortonEquivalentResistance;
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
}

void SimpleBranch::readNodeVoltage(TVectorD& nodeVoltageArray)
{
	for(int i=0;i<2;i++)
	{
		if(nodeNumber[i]==0)
			nodeVoltage[i]=0;
		else
			nodeVoltage[i]=nodeVoltageArray[nodeNumber[i]];
	}
}

void SimpleBranch::calculateBranchVoltage()
{
	branchVoltage_1 = branchVoltage;
	branchVoltage = nodeVoltage[0]-nodeVoltage[1];
}

void SimpleBranch::calculateBranchCurrent()
{
	branchCurrent_1 = branchCurrent;
	branchCurrent = nortonEquivalentCurrent + branchVoltage/nortonEquivalentResistance;
}

void SimpleBranch::formNodeNortonEquivalentCurrentArray(TVectorD& nodeNortonEquivalentCurrentArray)
{
	int from = nodeNumber[0];
	int to = nodeNumber[1];
	if(isSeriesOrNot)
	{
		nodeNortonEquivalentCurrentArray(from)-=nortonEquivalentCurrent;
		nodeNortonEquivalentCurrentArray(to)  +=nortonEquivalentCurrent;
	}	
	else
	{
		if(to==0) nodeNortonEquivalentCurrentArray(from)-=nortonEquivalentCurrent;
		else nodeNortonEquivalentCurrentArray(to)+=nortonEquivalentCurrent;
	}
}

void SimpleBranch::formConductanceMatrix(TMatrixD& conductanceMatrix)
{
	int from = nodeNumber[0];
	int to = nodeNumber[1];
	double nortonEquivalentConductor = 1/nortonEquivalentResistance;
	if(isSeriesOrNot)
	{
		conductanceMatrix(to,to)+=nortonEquivalentConductor;
		conductanceMatrix(from,from)+=nortonEquivalentConductor;
		conductanceMatrix(from,to)-=nortonEquivalentConductor;
		conductanceMatrix(to,from)-=nortonEquivalentConductor;
	}
	else if(from==0)
		conductanceMatrix(to,to)+=nortonEquivalentConductor;
	else
		conductanceMatrix(from,from)+=nortonEquivalentConductor;
}


void SimpleBranch::saveBranchCurrent(TMatrixD& branchCurrentMatrix,int& ptr,int counter)
{
	branchCurrentMatrix(counter,ptr) = branchCurrent;
	ptr++;
}

void SimpleBranch::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	branchCurrentMatrix_1[counter-1][ptr] = branchCurrent;
	ptr++;
}

//插值函数，普通元件只对电压和电流进行插值
void SimpleBranch::interpolate(double ratio)
{
	branchCurrent = (1-ratio)*branchCurrent_1 + ratio*branchCurrent;
	branchVoltage = (1-ratio)*branchVoltage_1 + ratio*branchVoltage;
	nortonEquivalentCurrent = (1-ratio)*nortonEquivalentCurrent_1 + ratio*nortonEquivalentCurrent;
}

//更新开关处理过程中用于存储结果的变量
void SimpleBranch::updateResult(int updateTypeNum)
{
	double V_temp,I_temp;//临时变量，交换数据时使用

	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		branchCurrent_2 = branchCurrent_1;
		branchVoltage_2 = branchVoltage_1;
		nortonEquivalentCurrent_2=nortonEquivalentCurrent_1;
		break;
	case 2://将_2变量的数值存入_1变量中
		branchCurrent_1 = branchCurrent_2;
		branchVoltage_1 = branchVoltage_2;
		nortonEquivalentCurrent_1=nortonEquivalentCurrent_2;
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

// 配合平均模型的迭代校正算法使用
// 目前仅PWMConverter和SimpleBranch中添加了这两个函数，xuyin,20121208
// 存储内部变量，以便在迭代校正时恢复
void SimpleBranch::storeInternalVariables()
{
	branchVoltage_bak = branchVoltage;
	branchCurrent_bak = branchCurrent;
	nortonEquivalentCurrent_bak = nortonEquivalentCurrent;
}
// 迭代校正时恢复内部变量
void SimpleBranch::restoreInternalVariables()
{
	branchVoltage = branchVoltage_bak;
	branchCurrent = branchCurrent_bak;
	nortonEquivalentCurrent = nortonEquivalentCurrent_bak;
}