#include "TimeSwitch.h"

#include <iostream>
using namespace std;

TimeSwitch::TimeSwitch(int id,int fromNode,int toNode,double onValue,double offValue,int state,int actionTimeCase)
{
	type = 8;
	this->id = id;
	nodeNumber[0] = fromNode;
	nodeNumber[1] = toNode;
	this->onValue = onValue;
	this->offValue = offValue;
	this->state = state;
	isSeriesOrNot=0;
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;

	actionTimeArray.ResizeTo(1,2);
	switch(actionTimeCase)
	{
		
	case 1:
		actionTimeArray(1)=3.0;
		actionTimeArray(2)=4.0;
		break;
	case 2:
		actionTimeArray(1)=1020e-3;
		break;
	default:
		cout<<"Warning：actionTimeCase"<<actionTimeCase<<"does not exist!"<<endl;
	}
	
}


void TimeSwitch::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
}

void TimeSwitch::calculateNortonEquivalentResistance(double time)
{//计算支路的诺顿等效电阻
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;
}

bool TimeSwitch::checkSwitch(double time)
{
	for(int i=1;i<=actionTimeArray.GetNrows();i++)
	{
		if(actionTimeArray(i)>time-deltaT/2 && actionTimeArray(i)<time+deltaT/2)
		{//actionTimeArray(i)==time
			if(state == 1){state = 0;calculateNortonEquivalentResistance(time);}
			else {state = 1;calculateNortonEquivalentResistance(time);}
			return 1;
		}
	}
		return 0;
}
