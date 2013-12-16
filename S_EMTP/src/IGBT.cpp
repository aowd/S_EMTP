#include "IGBT.h"
#include <cmath>
#include <iostream>
using namespace std;


IGBT::IGBT(int id,int fromNode,int toNode,double onValue,double offValue,int state,int actionTimeCase)
{
	type = 13;
	this->id = id;
	nodeNumber[0] = fromNode;
	nodeNumber[1] = toNode;
	this->onValue = onValue;
	this->offValue = offValue;
	this->state = state;
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;

	IGBT_number=actionTimeCase;

}

void IGBT::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent_1 = 0;
	nortonEquivalentCurrent = 0;
}

void IGBT::calculateNortonEquivalentResistance(double time)
{//计算支路的诺顿等效电阻
	if(state==0)
		nortonEquivalentResistance = offValue;
	else
		nortonEquivalentResistance = onValue;
}

bool IGBT::checkSwitch(double time)
{//检测开关动作
	int state_last;
	state_last=state;
	int tmp_state;
	state=switchController(time,IGBT_number);//更新当前状态
	//状态发生该变，则重新形成导纳阵
	if (state != state_last )
	{
		calculateNortonEquivalentResistance(time);
		return 1;
	}
	else
	{
		return 0;
	}

}

bool IGBT::switchController(double time,int SwitchNnumber)
{
	double time_last;
	double time_next;
	int sate_last;
	int sate_next;
	int tmp_state;

	time_next=time+5.0e-7;
	sate_next=camparePWM(time_next,SwitchNnumber);
	tmp_state=camparePWM(time,SwitchNnumber);
	////////////////////////////////==========//////////////////////////////////////////
	/////如果当前的状态和很短时间之后的状态不同，则说明状态马上要变了，那就提前变/////
	////////////////////////////////==========//////////////////////////////////////////
	if (tmp_state!=sate_next)
	{
		tmp_state=sate_next;
	}

	/////step2:根据1号管子的状态判断其他管子的状态
	if (SwitchNnumber%2==0)
	{
		tmp_state=1-tmp_state;
	}
	
	return tmp_state;
}

bool IGBT::camparePWM(double time,int SwitchNnumber)
{
	int state_1;
	double frequencyCW;
	double VmaxCW;
	double VoltageCW;
	double frequencyRef;
	double VmaxRef;
	double thet[60]={0,0,180,180,-72,-72,108,108,-144,-144,36,36,-216,-216,-36,-36,-288,-288,-108,-108,-12,-12,168,168,-84,-84,96,96,-156,-156,24,24,-228,-228,-48,-48,-300,-300,-120,-120,-24,-24,156,156,-96,-96,84,84,-168,-168,12,12,-240,-240,-60,-60,-312,-312,-132,-132};
	double VoltageRef;

	frequencyCW=2000.0;
	VmaxCW=2500.0;
	VmaxRef=2500.0*sqrt(2.0);
	frequencyRef=20.0;

	VoltageRef=VmaxRef*sin(2.0*PI*frequencyRef*time+thet[SwitchNnumber-1]*PI/180);

	double t=time;
	while(t>(1.0/frequencyCW))
	{
		t=t-1.0/frequencyCW;		
	}

	if (t<=0.5/frequencyCW)
	{
		/*VoltageCW = 2.0*frequencyCW*VmaxCW*t;*/         //for VoltageCW:0~VmaxCW
		VoltageCW = (4.0*frequencyCW*t-1)*VmaxCW;         //for VoltageCW:-VmaxCW=0~VmaxCW

	}
	else if (t>0.5/frequencyCW)
	{
		/*VoltageCW =(-2.0*frequencyCW*t+2.0)*VmaxCW;*/   //for VoltageCW:0~VmaxCW
		VoltageCW = (-4.0*frequencyCW*t+3.0)*VmaxCW;      //for VoltageCW:-VmaxCW=0~VmaxCW
	}

	if (VoltageRef>VoltageCW)
	{
		state_1=1;
	}
	else
	{
		state_1=0;
	}

	return state_1;
}
