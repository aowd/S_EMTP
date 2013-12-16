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

//��ֵ��������ͨԪ��ֻ�Ե�ѹ�͵������в�ֵ
void SimpleBranch::interpolate(double ratio)
{
	branchCurrent = (1-ratio)*branchCurrent_1 + ratio*branchCurrent;
	branchVoltage = (1-ratio)*branchVoltage_1 + ratio*branchVoltage;
	nortonEquivalentCurrent = (1-ratio)*nortonEquivalentCurrent_1 + ratio*nortonEquivalentCurrent;
}

//���¿��ش�����������ڴ洢����ı���
void SimpleBranch::updateResult(int updateTypeNum)
{
	double V_temp,I_temp;//��ʱ��������������ʱʹ��

	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		branchCurrent_2 = branchCurrent_1;
		branchVoltage_2 = branchVoltage_1;
		nortonEquivalentCurrent_2=nortonEquivalentCurrent_1;
		break;
	case 2://��_2��������ֵ����_1������
		branchCurrent_1 = branchCurrent_2;
		branchVoltage_1 = branchVoltage_2;
		nortonEquivalentCurrent_1=nortonEquivalentCurrent_2;
		break;
	default:
		cerr<<"��������ȷ�ĸ������ͱ�ţ�"<<endl;
		exit(1);
	}
}

// ���ƽ��ģ�͵ĵ���У���㷨ʹ��
// Ŀǰ��PWMConverter��SimpleBranch�������������������xuyin,20121208
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
void SimpleBranch::storeInternalVariables()
{
	branchVoltage_bak = branchVoltage;
	branchCurrent_bak = branchCurrent;
	nortonEquivalentCurrent_bak = nortonEquivalentCurrent;
}
// ����У��ʱ�ָ��ڲ�����
void SimpleBranch::restoreInternalVariables()
{
	branchVoltage = branchVoltage_bak;
	branchCurrent = branchCurrent_bak;
	nortonEquivalentCurrent = nortonEquivalentCurrent_bak;
}