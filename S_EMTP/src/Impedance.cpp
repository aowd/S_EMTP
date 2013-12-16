#include "Impedance.h"
#include <iostream>
using namespace std;

Impedance::Impedance(int id,int fromNode,int toNode,double resistance,double inductance)
{
	type=21;
	this->id=id;
	nPort=2;
	nodeNumber=new int[2];
	nodeNumber[0]=fromNode;
	nodeNumber[1]=toNode;
	branchCurrent=0;
	branchCurrent_1=0;
	branchCurrent_2=0;
	branchVoltage=0;
	isUserDef=0;
	need_NEC=1;
	this->resistance=resistance;
	this->inductance=inductance;
	nodeVoltage=new double[2];
	nodeVoltage[0]=0;
	nodeVoltage[1]=0;
	nortonEquivalentCurrent=0;
	nortonEquivalentCurrent_1=0;
	nortonEquivalentCurrent_2=0;
	nortonEquivalentResistance=resistance+2*inductance/deltaT;
	if (fromNode==0&&toNode==0)
	{
		cerr<<"֧·���˵㲻�ܾ��ӵ�"<<endl;
	}
}
Impedance::~Impedance()
{
	delete[] nodeNumber,nodeVoltage;
}
void Impedance::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)//��ʼ��֧·��ѹ����
{
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();
	branchCurrent=initialCurrentArray(ptr);
	ptr++;
}
void Impedance::readNodeVoltage(TVectorD& nodeVoltageArray)//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
{
	for (int i=0;i<2;i++)
	{
		if (nodeNumber[i]==0)
		{
			nodeVoltage[i]=0;
		}
		else
		{
			nodeVoltage[i]=nodeVoltageArray(nodeNumber[i]);
		}
	}
}
void Impedance::calculateBranchVoltage()//����֧·��ѹ
{
	branchVoltage_1 = branchVoltage;
	branchVoltage=nodeVoltage[0]-nodeVoltage[1];
}
void Impedance::calculateBranchCurrent()//����֧·����
{
	branchCurrent_1=branchCurrent;
	branchCurrent=branchVoltage/nortonEquivalentResistance+nortonEquivalentCurrent;
}
void Impedance::calculateNortonEquivalentCurrent(double time)//����֧·��ŵ�ٵ�Ч��·�еĵ�����
{
	nortonEquivalentCurrent_1=nortonEquivalentCurrent;
	nortonEquivalentCurrent=(branchVoltage+(2*inductance/deltaT-resistance)*branchCurrent)/nortonEquivalentResistance;
}
void Impedance::calculateNortonEquivalentResistance(double time)//����֧·��ŵ�ٵ�Ч����
{
	nortonEquivalentResistance=resistance+2*inductance/deltaT;
}
void Impedance::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)//�γɽڵ�ŵ�ٵ�Ч��������
{
	for (int i=0;i<2;i++)
	{
		if (nodeNumber[i]!=0)
		{
			if (i==0)
			{
				nodeNortonEquivalentCurrentArray(nodeNumber[0])-=nortonEquivalentCurrent;
			}
			else
				nodeNortonEquivalentCurrentArray(nodeNumber[1])+=nortonEquivalentCurrent;
		} 
	}
}
void Impedance::formConductanceMatrix(TMatrixD &conductanceMatrix)//�γɽڵ㵼����
{
	int from=nodeNumber[0];
	int to=nodeNumber[1];
	double nortonEquivalentConductance=1/nortonEquivalentResistance;
	if (from!=0&&to!=0)
	{
		conductanceMatrix(from,from)+=nortonEquivalentConductance;
		conductanceMatrix(to,to)+=nortonEquivalentConductance;
		conductanceMatrix(from,to)-=nortonEquivalentConductance;
		conductanceMatrix(to,from)-=nortonEquivalentConductance;
	} 
	else
	{
		if (from==0)
		{
			conductanceMatrix(to,to)+=nortonEquivalentConductance;
		} 
		else
		{
			conductanceMatrix(from,from)+=nortonEquivalentConductance;
		}
	}
}
void Impedance::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)//����֧·����
{
	branchCurrentMatrix(counter,ptr)=branchCurrent;
	ptr++;
}
void Impedance::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)//����֧·����
{
	branchCurrentMatrix_1[counter-1][ptr]=branchCurrent;
	ptr++;
}
double Impedance::getNortonEquivalentResistance(){return nortonEquivalentResistance;}
double Impedance::getNortonEquivalentCurrent(){return nortonEquivalentCurrent;}
double Impedance::getNortonEquivalentCurrent_1(){return nortonEquivalentCurrent_1;}
double Impedance::getNortonEquivalentCurrent_2(){return nortonEquivalentCurrent_2;}