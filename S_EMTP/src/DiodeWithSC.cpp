#include "DiodeWithSC.h"
#include <iostream>
using namespace std;

DiodeWithSC::DiodeWithSC(int id,int fromNode,int toNode,double onValue,double offValue,double Rd,double Cd,int state)
{
	type = 77;
	this->id = id;
	nodeNumber[0] = fromNode;
	nodeNumber[1] = toNode;
	this->onValue = onValue;
	this->offValue = offValue;
	this->Rd = Rd;
	this->Cd = Cd;
	this->state = state;
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;
	forwardVoltageDrop = 0;
	forwardBreakoverVoltage = 1e8;
	reverseWithstandVoltage = 1e8;
}

void DiodeWithSC::calculateNortonEquivalentCurrent(double time)
{
	double iC,vC,iC_norton,Rv;
	
	if (state == 0)
	{
		Rv = offValue;
	}else
	{
		Rv = onValue;
	}

	iC = branchCurrent - branchVoltage/Rv;
	vC = branchVoltage - iC*Rd;
	iC_norton = iC + (2*Cd/deltaT)*vC;
	
	nortonEquivalentCurrent_1 = nortonEquivalentCurrent;
	nortonEquivalentCurrent = (1/Rd)/(1/Rd+2*Cd/deltaT)*iC_norton;
}

void DiodeWithSC::calculateNortonEquivalentResistance(double time)
{
	if(state==0)
		nortonEquivalentResistance = offValue*(Rd+deltaT/2/Cd)/(offValue+Rd+deltaT/2/Cd);
	else
		nortonEquivalentResistance = onValue*(Rd+deltaT/2/Cd)/(onValue+Rd+deltaT/2/Cd);
}

bool DiodeWithSC::checkSwitch(double time)
{
	if(state==0 && branchVoltage>forwardVoltageDrop)
	{
		state = 1;
		calculateNortonEquivalentResistance(time);
		return 1;
	}
	else if(state==1 && branchVoltage<0)
	{
		state = 0;
		calculateNortonEquivalentResistance(time);
		return 1;
	}
	else
		return 0;
}