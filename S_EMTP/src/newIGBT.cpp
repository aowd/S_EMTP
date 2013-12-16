#include "newIGBT.h"
#include <iostream>
using namespace std;

newIGBT::newIGBT(int id,int fromNode,int toNode,double onValue,double offValue,int CtrlSystemNode)
{
	type = 22;
	this->id = id;
	nodeNumber[0] = fromNode;
	nodeNumber[1] = toNode;
	this->onValue = onValue;
	this->offValue = offValue;
	this->CtrlSystemNode=CtrlSystemNode;
	forwardVoltageDrop = 1e-5;
	//forwardVoltageDrop = 0;
	forwardBreakoverVoltage = 1e8;
	reverseWithstandVoltage = 1e8;
	state=0;
	nortonEquivalentResistance=(state==1)?onValue:offValue;
	controlChangeOrNot=0;
	isSeriesOrNot = isSeries();
	controlSignal=0;
	controlSignal_1=controlSignal;
	controlSignal_2=controlSignal_1;
	switchMode = 0;
}

void newIGBT::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
}

void newIGBT::calculateNortonEquivalentResistance(double time)
{
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;
}

void newIGBT::switchIt()
{
	if (state == 1)
	{
		state = 0;
		//cout<<"IGBT"<<id<<"关断"<<endl;
	} 
	else
	{
		state = 1;
		//cout<<"IGBT"<<id<<"导通"<<endl;
	}
	calculateNortonEquivalentResistance(1);
}

void newIGBT::modifyConductanceMatrix(TMatrixD &conductanceMatrix)
{
	double temp;//代表诺顿等值导纳的增量
	if (state == 1)
	{
		temp = 1/onValue-1/offValue;
	}
	else
	{
		temp = 1/offValue-1/onValue;
	}

	int from = nodeNumber[0];
	int to = nodeNumber[1];
	if(isSeriesOrNot)
	{
		conductanceMatrix(to,to)+=temp;
		conductanceMatrix(from,from)+=temp;
		conductanceMatrix(from,to)-=temp;
		conductanceMatrix(to,from)-=temp;
	}
	else if(from==0)
		conductanceMatrix(to,to)+=temp;
	else
		conductanceMatrix(from,from)+=temp;
}

//判断开关是否需要动作
bool newIGBT::checkSwitch(double time)
{
	controlChangeOrNot=!(controlSignal_2==controlSignal_1);
	if ((controlSignal_1==1 && branchVoltage>forwardVoltageDrop && state==0)||((controlSignal_1==0 ||  branchCurrent<-1e-5)&& state==1))
	//if ((controlSignal_1==1 && branchVoltage>forwardVoltageDrop && state==0)||((controlSignal_1==0 ||  branchCurrent<0)&& state==1))
	{
		return 1;
	}
	//controlChangeOrNot=!(controlSignal==controlSignal_1);
	//if ((controlSignal==1 && branchVoltage>forwardVoltageDrop && state==0)||((controlSignal==0 ||  branchCurrent<0)&& state==1))
	//{
	//	return 1;
	//}
	else
	{
		return 0;
	}
}

//返回开关动作点的插值比
double newIGBT::getSwitchRatio()
{
	if ((controlChangeOrNot&&branchVoltage_1>0&&state==0)||(controlChangeOrNot&&branchCurrent>0&&state==1))
	{
		switchMode = 1;
		return 0;
	} 
	//if ((controlChangeOrNot&&branchVoltage_1>0&&state==0)||(controlChangeOrNot&&branchCurrent>0&&state==1))
	//{
	//	switchMode = 1;
	//	return 1;
	//} 
	else	
	{
		switchMode = 0;
		if (state == 1)
		{
			return branchCurrent_1/(branchCurrent_1-branchCurrent);
		} 
	else
		{
			return branchVoltage_1/(branchVoltage_1-branchVoltage);
		}
	}
}

void newIGBT::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();
	branchCurrent = initialCurrentArray[ptr];
	ptr++;

	if (branchVoltage/branchCurrent<2*onValue)
	{
		state = 1;
		controlSignal=1;
		calculateNortonEquivalentResistance(time);
	}
	else
	{
		state = 0;
		controlSignal=0;
		calculateNortonEquivalentResistance(time);
	}
}

int newIGBT::getSwitchMode(){
	return switchMode;
}

void newIGBT::setControlledVariable(double* ctrlNodeValue)
{
	controlSignal_2=controlSignal_1;
	controlSignal_1=controlSignal;

	//controlSignal = (ctrlNodeValue[CtrlSystemNode-1]>0)?1:0;
	controlSignal = ctrlNodeValue[CtrlSystemNode-1];
}

void newIGBT::setControlledVariableForSwitch(double* ctrlNodeValue)
{
	//controlSignal = (ctrlNodeValue[CtrlSystemNode-1]>0)?1:0;
	controlSignal = ctrlNodeValue[CtrlSystemNode-1];
}
