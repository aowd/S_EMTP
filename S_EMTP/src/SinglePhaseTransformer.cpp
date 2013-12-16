#include "SinglePhaseTransformer.h"

#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

SinglePhaseTransformer::SinglePhaseTransformer(int id,int in_nodeNumber[],double apparentPowerRating,double voltageRating[]){//构造函数
	type = 5;
	nPort = 4;
	need_NEC=1;
	nodeNumber=new int[nPort];
	for(int i=0;i<nPort;i++){
		this->nodeNumber[i] = in_nodeNumber[i];
		this->nodeVoltageArray[i] = 0;
	}
	this->apparentPowerRating = apparentPowerRating;
	for(int i=0;i<2;i++)
		this->voltageRating[i] = voltageRating[i];
	frequency = 50;
	leakageReactance = 0.1;
	magnetizingCurrent = 1;
	coreAspectRatio[0] = 2;
	coreAspectRatio[1] = 1;

	for(int i=0;i<2;i++){
		branchCurrentArray[i] = 0;
		branchVoltageArray[i] = 0;
		fluxArray[i] = 0;
		nortonEquivalentCurrent[i] = 0;
		branchVoltageArray_1[i] = 0;
		branchCurrentArray_1[i] = 0;
	}

	calculatePermeanceMatrix();
	calculateAdmittanceMatrix();

	double fm = (permeance[0]*permeance[0]-permeance[1]*permeance[1])*voltageRating[0]*voltageRating[1];
	invert_MN[0] = permeance[0]*voltageRating[1]/fm;
	invert_MN[1] = -permeance[1]*voltageRating[1]/fm;
	invert_MN[2] = -permeance[1]*voltageRating[0]/fm;
	invert_MN[3] = permeance[0]*voltageRating[0]/fm;
}

void SinglePhaseTransformer::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
//初始化支路电压电流
{
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	branchCurrentArray[0] = initialCurrentArray(ptr+1);
	branchCurrentArray[1] = initialCurrentArray(ptr+2);
	ptr+=2;
}
void SinglePhaseTransformer::readNodeVoltage(TVectorD& nodeVoltageArray)//从节点电压数组中读支路两节点电压
{
	for(int i=0;i<4;i++)
	{
		if(nodeNumber[i]==0)
			this->nodeVoltageArray[i]=0;
		else
			this->nodeVoltageArray[i]=nodeVoltageArray[nodeNumber[i]];
	}
}
void SinglePhaseTransformer::calculateBranchVoltage()//计算支路电压
{
	branchVoltageArray_1[0] = branchVoltageArray[0];
	branchVoltageArray_1[1] = branchVoltageArray[1];
	branchVoltageArray[0] = nodeVoltageArray[0]-nodeVoltageArray[1];
	branchVoltageArray[1] = nodeVoltageArray[2]-nodeVoltageArray[3];

	fluxArray_1[0] = fluxArray[0];
	fluxArray_1[1] = fluxArray[1];
	fluxArray[0] = fluxArray[0] + (deltaT/2)*(branchVoltageArray_1[0]+branchVoltageArray[0])/voltageRating[0];
	fluxArray[1] = fluxArray[1] + (deltaT/2)*(branchVoltageArray_1[1]+branchVoltageArray[1])/voltageRating[1];

}

void SinglePhaseTransformer::calculateBranchCurrent()//计算支路电流
{
	branchCurrentArray_1[0] = branchCurrentArray[0];
	branchCurrentArray_1[1] = branchCurrentArray[1];

	branchCurrentArray[0] = conductance[0]*branchVoltageArray[0]
		+ conductance[1]*branchVoltageArray[1] + nortonEquivalentCurrent[0];
	branchCurrentArray[1] = conductance[1]*branchVoltageArray[0]
		+ conductance[2]*branchVoltageArray[1] + nortonEquivalentCurrent[1];
	//cout<<branchCurrentArray_1[0]<<endl;//for test 09.05.19
	//cout<<branchCurrentArray_1[1]<<endl;//for test 09.05.19
}

void SinglePhaseTransformer::calculateNortonEquivalentCurrent(double time)//计算支路的诺顿等效电路中的电流项
{
	nortonEquivalentCurrent_1[0] = nortonEquivalentCurrent[0];
	nortonEquivalentCurrent_1[1] = nortonEquivalentCurrent[1];
	nortonEquivalentCurrent[0] = conductance[0]*branchVoltageArray[0] + conductance[1]*branchVoltageArray[1]
		+ invert_MN[0]*fluxArray[0] + invert_MN[1]*fluxArray[1];
	nortonEquivalentCurrent[1] = conductance[1]*branchVoltageArray[0] + conductance[2]*branchVoltageArray[1]
		+ invert_MN[2]*fluxArray[0] + invert_MN[3]*fluxArray[1];
		//cout<<nortonEquivalentCurrent[0]<<nortonEquivalentCurrent[1]<<endl;
}

void SinglePhaseTransformer::calculatePermeanceMatrix()//计算铁心磁导矩阵
{
	double p12 = 200*(1+coreAspectRatio[0]/coreAspectRatio[1])/(2*PI*frequency*apparentPowerRating*magnetizingCurrent);
	double p3 = 100*(1+coreAspectRatio[1]/coreAspectRatio[0])/(2*PI*frequency*apparentPowerRating*magnetizingCurrent);
	double p45 = leakageReactance/(4*PI*frequency*apparentPowerRating);

	double fz1 = ((2*p3+p45)*p45+p3*p12+p45*p12)*p12;
	double fz2 = p3*p12*p12;
	double fm = p12*p12+2*p12*p45+2*p3*p12+2*p3*p45+p45*p45;

	permeance[0] = fz1/fm;
	permeance[1] = fz2/fm;
}

void SinglePhaseTransformer::calculateAdmittanceMatrix()//计算诺顿等效导纳矩阵
{
	double N1 = voltageRating[0];
	double N2 = voltageRating[1];
	double M1 = permeance[0];
	double M2 = permeance[1];

	double c = deltaT/(2*(M1*M1-M2*M2));

	conductance[0] = c*M1/(N1*N1);
	conductance[1] = -c*M2/(N1*N2);
	conductance[2] = c*M1/(N2*N2);
}

void SinglePhaseTransformer::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
//形成节点诺顿等效电流向量
{
	int N1 = nodeNumber[0];
	int N2 = nodeNumber[2];

	nodeNortonEquivalentCurrentArray(N1) -= nortonEquivalentCurrent[0];
	nodeNortonEquivalentCurrentArray(N2) -= nortonEquivalentCurrent[1];
}

void SinglePhaseTransformer::formConductanceMatrix(TMatrixD& conductanceMatrix)//形成节点导纳阵
{
	int N1 = nodeNumber[0];
	int N2 = nodeNumber[2];
	//cout<<N1<<"\t"<<N2<<endl;
	conductanceMatrix(N1,N1) += conductance[0];
	conductanceMatrix(N1,N2) += conductance[1];
	conductanceMatrix(N2,N1) += conductance[1];
	conductanceMatrix(N2,N2) += conductance[2];
	//conductanceMatrix.Print();

}

void SinglePhaseTransformer::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)//保存支路电流
{
	branchCurrentMatrix(counter,ptr) = branchCurrentArray[0];
	branchCurrentMatrix(counter,ptr+1) = branchCurrentArray[1];
	ptr+=2;
}

void SinglePhaseTransformer::interpolate(double ratio)//给定比值，对支路的电压电流进行插值
{
	branchCurrentArray[0] = (1-ratio)*branchCurrentArray_1[0] + ratio*branchCurrentArray[0];
	branchCurrentArray[1] = (1-ratio)*branchCurrentArray_1[1] + ratio*branchCurrentArray[1];
	branchVoltageArray[0] = (1-ratio)*branchVoltageArray_1[0] + ratio*branchVoltageArray[0];
	branchVoltageArray[1] = (1-ratio)*branchVoltageArray_1[1] + ratio*branchVoltageArray[1];
	nortonEquivalentCurrent[0]= (1-ratio)*nortonEquivalentCurrent_1[0]+ ratio*nortonEquivalentCurrent[0];
    nortonEquivalentCurrent[1]= (1-ratio)*nortonEquivalentCurrent_1[1]+ ratio*nortonEquivalentCurrent[1];
	fluxArray[0] = (1-ratio)*fluxArray_1[0]+ratio*fluxArray[0];
	fluxArray[1] = (1-ratio)*fluxArray_1[1]+ratio*fluxArray[1];
}