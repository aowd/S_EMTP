#include "Inverter.h"
#include <cmath>
#include <iostream>
using namespace std;

Inverter::Inverter(int id,int DCNode_first,int DCNode_last,int ACNode_first,int ACNode_last,int Inverter_number,double onValue,double offValue)
{
	type = 11;
	this->id = id;
	nPort=12;
	nodeNumber=new int[nPort];

	this->onValue = onValue;
	this->offValue = offValue;

	for (int i=0;i<10;i++)
	{
		nodeNumber[i]=ACNode_first+i;
	}
	nodeNumber[10] = DCNode_first;
	nodeNumber[11] = DCNode_last;


	for(int i=0;i<10;i++)
	{
		this->IGBTstate[2*i] = 0;
		this->IGBTstate[2*i+1] = 0;
	}

	for(int i=0;i<20;i++)
	{
		//this->IGBTstate[i] = 0;
		if(IGBTstate[i]==0)
			nortonEquivalentResistance[i] = offValue;
		else
			nortonEquivalentResistance[i] = onValue;
	}
	this->Inverter_number=Inverter_number;

	for (int k=0;k<20;k++)
	{
		IGBTVoltage[k] = 0;
		IGBTCurrent[k] = 0;
		IGBTCurrent_1[k] = 0;
		IGBTCurrent_2[k] = 0;
	}
//
//#ifdef  WIN32
//	ifstream infile("configure.txt",ios::in);
//#else
//	ifstream infile("../configure.txt",ios::in);
//#endif
//	double tmp_finishTime,tmp_deltaT;
//	infile>>tmp_finishTime>>tmp_deltaT;
//	infile.close();   
//	int nstep=tmp_finishTime/tmp_deltaT+10;//留一点裕量出来
//	for (int i=0;i<10;i++)
//	{
//		inverterState[i]=new bool[nstep];
//	}
//
//	double deltaT = 50e-6;
//	for (int k=0;k<10;k++)
//	{
//		SwitchNnumber=20*(Inverter_number-1)+2*k+1;
//		for (int i=1;i<=nstep;i++)
//		{
//			inverterState[k][i-1]=switchController(deltaT*i,SwitchNnumber);
//		}
//	}
}

void Inverter::calculateNortonEquivalentCurrent(double time)
{
	nortonEquivalentCurrent = 0;
}

void Inverter::calculateNortonEquivalentResistance(double time)
{//计算支路的诺顿等效电阻
	for (int i=0;i<20;i++)
	{
		if(IGBTstate[i]==0)
			nortonEquivalentResistance[i] = offValue;
		else
			nortonEquivalentResistance[i] = onValue;
	}
}

bool Inverter::checkSwitch(int counter,TMatrixD &conductanceMatrix)
{//检测开关动作
	int state_last;
	int tmp_state;
	int statechange=0;
	int SwitchNnumber;
	
	for (int i=0;i<10;i++)
	{
		SwitchNnumber=20*(Inverter_number-1)+2*i+1;
		state_last=IGBTstate[2*i];
		SwitchNnumber=20*(Inverter_number-1)+2*i+1;
		IGBTstate[2*i]=switchController(deltaT*(counter-1),SwitchNnumber);//更新当前状态
		IGBTstate[2*i+1]=1-IGBTstate[2*i];
		//状态发生该变，则重新形成导纳阵
		if (counter>2 && IGBTstate[2*i] != state_last )
		{
			statechange++;
			modifyConductanceMatrix(conductanceMatrix,i);		
		}
		else if (counter==2)
		{
			statechange++;
			modifyConductanceMatrixforfirststep(conductanceMatrix,i);
		}
	}
			
	if (statechange==0)
	{
		return 0;
	}
	else
	{
		calculateNortonEquivalentResistance(1);
		return 1;
	}
}

void Inverter::modifyConductanceMatrix(TMatrixD &conductanceMatrix,int k)
{
	int k1,k2,k3;
	k1 = nodeNumber[k];
	k2 = nodeNumber[10];
	k3 = nodeNumber[11];
	
	double temp;
	if (IGBTstate[2*k]==1)
	{
		temp = 1/onValue-1/offValue;
	}
	else
	{
		temp = 1/offValue-1/onValue;
	}

	conductanceMatrix(k2,k2) += temp;
	conductanceMatrix(k2,k1) -= temp;
	conductanceMatrix(k1,k2) -= temp;

	conductanceMatrix(k3,k3) -= temp;
	conductanceMatrix(k3,k1) += temp;
	conductanceMatrix(k1,k3) += temp;
}

void Inverter::modifyConductanceMatrixforfirststep(TMatrixD &conductanceMatrix,int k)
{
	int k1,k2,k3;
	k1 = nodeNumber[k];
	k2 = nodeNumber[10];
	k3 = nodeNumber[11];

	double temp;
	temp = 1/onValue-1/offValue;
	if (IGBTstate[2*k]==1)
	{
		conductanceMatrix(k2,k2) += temp;
		conductanceMatrix(k1,k1) += temp;
		conductanceMatrix(k2,k1) -= temp;
		conductanceMatrix(k1,k2) -= temp;
	}
	else
	{
		conductanceMatrix(k3,k3) += temp;
		conductanceMatrix(k1,k1) += temp;
		conductanceMatrix(k3,k1) -= temp;
		conductanceMatrix(k1,k3) -= temp;		
	}
}

bool Inverter::switchController(double time,int SwitchNnumber)
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
	//if (SwitchNnumber%2==0)
	//{
	//	tmp_state=1-tmp_state;
	//}
	return tmp_state;
}

bool Inverter::camparePWM(double time,int SwitchNnumber)
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
		VoltageCW = (4.0*frequencyCW*t-1)*VmaxCW;         //for VoltageCW:-VmaxCW~VmaxCW
	}
	else if (t>0.5/frequencyCW)
	{
		/*VoltageCW =(-2.0*frequencyCW*t+2.0)*VmaxCW;*/   //for VoltageCW:0~VmaxCW
		VoltageCW = (-4.0*frequencyCW*t+3.0)*VmaxCW;      //for VoltageCW:-VmaxCW~VmaxCW
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
void Inverter::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();
	for (int i=0;i<20;i++)
	{
		IGBTCurrent[i] = initialCurrentArray[ptr+i];
	}
	ptr+=20;
}

void Inverter::readNodeVoltage(TVectorD& nodeVoltageArray)
{
	for(int i=0;i<12;i++)
	{
		if(nodeNumber[i]==0)
			nodeVoltage[i]=0;
		else
			nodeVoltage[i]=nodeVoltageArray[nodeNumber[i]];
	}
}

void Inverter::calculateBranchVoltage()
{
	calculateIGBTVoltage();
}

void Inverter::calculateIGBTVoltage()
{
	for (int i=0;i<10;i++)
	{
		IGBTVoltage[2*i]=nodeVoltage[10]-nodeVoltage[i];
		IGBTVoltage[2*i+1]=nodeVoltage[i]-nodeVoltage[11];
	}
}

void Inverter::calculateBranchCurrent()
{
	calculateIGBTCurrent();
}

void Inverter::calculateIGBTCurrent()
{
	for (int i=0;i<20;i++)
	{
		IGBTCurrent_1[i] = IGBTCurrent[i];
		IGBTCurrent[i]=IGBTVoltage[i]/nortonEquivalentResistance[i];
	}
}

void Inverter::formNodeNortonEquivalentCurrentArray(TVectorD& nodeNortonEquivalentCurrentArray)
{
}

void Inverter::formConductanceMatrix(TMatrixD& conductanceMatrix)
{
	int from,to;
	double nortonEquivalentConductor;
	for (int i=0;i<10;i++)
	{
		from=nodeNumber[10];
		to=nodeNumber[i]; 
		nortonEquivalentConductor = 1/nortonEquivalentResistance[2*i];
		conductanceMatrix(to,to)+=nortonEquivalentConductor;
		conductanceMatrix(from,from)+=nortonEquivalentConductor;
		conductanceMatrix(from,to)-=nortonEquivalentConductor;
		conductanceMatrix(to,from)-=nortonEquivalentConductor;
		
		from=to;
		to=nodeNumber[11];
		nortonEquivalentConductor = 1/nortonEquivalentResistance[2*i+1];
		conductanceMatrix(to,to)+=nortonEquivalentConductor;
		conductanceMatrix(from,from)+=nortonEquivalentConductor;
		conductanceMatrix(from,to)-=nortonEquivalentConductor;
		conductanceMatrix(to,from)-=nortonEquivalentConductor;
	}
}


void Inverter::saveBranchCurrent(TMatrixD& branchCurrentMatrix,int& ptr,int counter)
{
	for (int i=0;i<20;i++)
	{
		branchCurrentMatrix(counter,ptr) = IGBTCurrent[i];
		ptr++;
	}
}

void Inverter::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int i=0;i<20;i++)
	{
		branchCurrentMatrix_1[counter-1][ptr] = IGBTCurrent[i];
		ptr++;
	}
}

//插值函数，普通元件只对电压和电流进行插值
void Inverter::interpolate(double ratio)
{
	for (int k=0;k<20;k++)
	{
		IGBTCurrent[k] = (1-ratio)*IGBTCurrent_1[k] + ratio*IGBTCurrent[k];
	}
}

//更新开关处理过程中用于存储结果的变量
void Inverter::updateResult(int updateTypeNum)
{
	double V_temp,I_temp;//临时变量，交换数据时使用

	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int k=0;k<20;k++)
		{
			IGBTCurrent_2[k] = IGBTCurrent_1[k];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int k=0;k<20;k++)
		{
			IGBTCurrent_1[k] = IGBTCurrent_2[k];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}