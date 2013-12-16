#include "ThreePhaseRectifierBridge.h"
#include "Component.h"
#include <iostream>
ThreePhaseRectifierBridge::ThreePhaseRectifierBridge(int id,int	*ACNode,int *DCNode,int *state,double onValue,double offValue)
{
	type = 17;
	nPort=5;
	nodeNumber=new int[nPort];
	this->id = id;
	nodeNumber[0] = ACNode[0];
	nodeNumber[1] = ACNode[1];
	nodeNumber[2] = ACNode[2];
	nodeNumber[3] = DCNode[0];
	nodeNumber[4] = DCNode[1];
	numberdiode();

	this->onValue = onValue;
	this->offValue = offValue;

	for (int i=0;i<6;i++)
	{
		diodeState[i]=state[i];//6个二极管的状态
	}
	//diodeState[4]=1;
	//diodeState[5]=1;
	bridgeMode=0;
	bridgeMode_tmp=4;
	
	forwardVoltageDrop = 0;
	forwardBreakoverVoltage = 1e8;
	reverseWithstandVoltage = 1e8;
}

void ThreePhaseRectifierBridge::numberdiode()
{
	diodenodeNo[0][0] = nodeNumber[0];
	diodenodeNo[0][1] = nodeNumber[3];
	diodenodeNo[1][0] = nodeNumber[4];
	diodenodeNo[1][1] = nodeNumber[2];
	diodenodeNo[2][0] = nodeNumber[1];
	diodenodeNo[2][1] = nodeNumber[3];
	diodenodeNo[3][0] = nodeNumber[4];
	diodenodeNo[3][1] = nodeNumber[0];
	diodenodeNo[4][0] = nodeNumber[2];
	diodenodeNo[4][1] = nodeNumber[3];
	diodenodeNo[5][0] = nodeNumber[4];
	diodenodeNo[5][1] = nodeNumber[1];
}


void ThreePhaseRectifierBridge::calculateBranchCurrent()   //????
{
	calculateDiodeCurrent();
	calculatebranchCurrentArray();
}

void ThreePhaseRectifierBridge::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	nortonEquivalentCurrent = 0;
}

void ThreePhaseRectifierBridge::calculateNortonEquivalentResistance(double time)
{//计算支路的诺顿等效电阻
	for(int i=0;i<6;i++)
	{
		if(diodeState[i]==0)
			nortonEquivalentResistance[i] = offValue;
		else
			nortonEquivalentResistance[i] = onValue;
	}
}

void ThreePhaseRectifierBridge::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{//初始化支路电压电流
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();


	for (int i=0;i<6;i++)
	{
		branchCurrentArray[i] = initialCurrentArray[ptr+i];
	}
	ptr+=6;
	nortonEquivalentCurrent=0;

}

void ThreePhaseRectifierBridge::readNodeVoltage(TVectorD& nodeVoltageArray)
{//从节点电压数组中读支路两节点电压
	for(int i=0;i<5;i++)
	{
		if(nodeNumber[i]==0)
			nodeVoltage[i]=0;
		else
			nodeVoltage[i]=nodeVoltageArray[nodeNumber[i]];
	}
}

void ThreePhaseRectifierBridge::calculateDiodeVoltage()//added by GongYuan  08.12.05
{
	//for (int i=0;i<6;i++)
	//{
	//	diodeVoltage_1[i] = diodeVoltage[i];
	//	//diodeVoltage[i] = nodeVoltage[diodenodeNo[i][0]-1]-nodeVoltage[diodenodeNo[i][1]-1];
	//}
	diodeVoltage[0] = nodeVoltage[0] - nodeVoltage[3];
	diodeVoltage[1] = nodeVoltage[4] - nodeVoltage[2];
	diodeVoltage[2] = nodeVoltage[1] - nodeVoltage[3];
	diodeVoltage[3] = nodeVoltage[4] - nodeVoltage[0];
	diodeVoltage[4] = nodeVoltage[2] - nodeVoltage[3];
	diodeVoltage[5] = nodeVoltage[4] - nodeVoltage[1];
}

void ThreePhaseRectifierBridge::calculateDiodeCurrent()//计算支路电流 added by GongYuan  08.12.05
{
	for (int i=0;i<6;i++)
	{
		//diodeCurrent_1[i] = diodeCurrent[i];
		diodeCurrent[i] = diodeVoltage[i]/nortonEquivalentResistance[i];
		//cout<<"diodeCurrent:"<<diodeCurrent[i]<<endl;
	}
}

void ThreePhaseRectifierBridge::calculateBranchVoltage()  //????
{//计算支路电压
	calculateDiodeVoltage();
}

void ThreePhaseRectifierBridge::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)//added by GongYuan  08.12.11
{//形成节点诺顿等效电流向量
	//for (int i=0;i<6;i++)
	//{
	//	int from = diodenodeNo[i][0];
	//	int to = diodenodeNo[i][1];
	//	nodeNortonEquivalentCurrentArray(from)-=nortonEquivalentCurrent;
	//	nodeNortonEquivalentCurrentArray(to)  +=nortonEquivalentCurrent;
	//}
}
void ThreePhaseRectifierBridge::formConductanceMatrix(TMatrixD &conductanceMatrix)//added by GongYuan  08.12.11
{//形成节点导纳阵
	double nortonEquivalentConductor[6];
	//for (int i=0;i<6;i++)
	//{
	//	int from = diodenodeNo[i][0];
	//	int to = diodenodeNo[i][1];
	//	nortonEquivalentConductor[i] = 1/nortonEquivalentResistance[i];
	//	conductanceMatrix(to,to)+=nortonEquivalentConductor[i];
	//	conductanceMatrix(from,from)+=nortonEquivalentConductor[i];
	//	conductanceMatrix(from,to)-=nortonEquivalentConductor[i];
	//	conductanceMatrix(to,from)-=nortonEquivalentConductor[i];
	//}

	for (int i=0;i<6;i++)
	{
		nortonEquivalentConductor[i] = 1.0/nortonEquivalentResistance[i]; 
	}
	conductanceMatrix(nodeNumber[0],nodeNumber[0])+=nortonEquivalentConductor[0]+nortonEquivalentConductor[3];
	conductanceMatrix(nodeNumber[1],nodeNumber[1])+=nortonEquivalentConductor[2]+nortonEquivalentConductor[5];
	conductanceMatrix(nodeNumber[2],nodeNumber[2])+=nortonEquivalentConductor[4]+nortonEquivalentConductor[1];
	conductanceMatrix(nodeNumber[3],nodeNumber[3])+=nortonEquivalentConductor[0]+nortonEquivalentConductor[2]+nortonEquivalentConductor[4];
	conductanceMatrix(nodeNumber[4],nodeNumber[4])+=nortonEquivalentConductor[1]+nortonEquivalentConductor[3]+nortonEquivalentConductor[5];
	conductanceMatrix(nodeNumber[0],nodeNumber[3])-=nortonEquivalentConductor[0];
	conductanceMatrix(nodeNumber[0],nodeNumber[4])-=nortonEquivalentConductor[3];
	conductanceMatrix(nodeNumber[1],nodeNumber[3])-=nortonEquivalentConductor[2];
	conductanceMatrix(nodeNumber[1],nodeNumber[4])-=nortonEquivalentConductor[5];
	conductanceMatrix(nodeNumber[2],nodeNumber[3])-=nortonEquivalentConductor[4];
	conductanceMatrix(nodeNumber[2],nodeNumber[4])-=nortonEquivalentConductor[1];
	conductanceMatrix(nodeNumber[3],nodeNumber[0])-=nortonEquivalentConductor[0];
	conductanceMatrix(nodeNumber[3],nodeNumber[1])-=nortonEquivalentConductor[2];
	conductanceMatrix(nodeNumber[3],nodeNumber[2])-=nortonEquivalentConductor[4];
	conductanceMatrix(nodeNumber[4],nodeNumber[0])-=nortonEquivalentConductor[3];
	conductanceMatrix(nodeNumber[4],nodeNumber[1])-=nortonEquivalentConductor[5];
	conductanceMatrix(nodeNumber[4],nodeNumber[2])-=nortonEquivalentConductor[1];

}
void ThreePhaseRectifierBridge::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//保存支路电流
	for (int i=0;i<6;i++)
	{
		branchCurrentMatrix(counter,ptr+i) = diodeCurrent[i];
	}
	ptr+=6;
}

bool ThreePhaseRectifierBridge::checkSwitch(double time)
{//检测开关动作
	double switchtime_1;
	//int bridgeMode_1;
	//int statechange[2]={0,0};
	statechange[0]=0;
	statechange[1]=0;
	bridgeMode_1=bridgeMode;
	//checkbridgemode();

	for (int i=0;i<6;i++)
	{
		if(diodeState[i]==0 && diodeVoltage[i]>forwardVoltageDrop)
		{
			diodeState[i] = 1;
			statechange[0]++;
		}
		else if(diodeState[i]==1 && diodeCurrent[i]<0)
		{
			diodeState[i] = 0;
			statechange[1]++;
		}
	}
	//checkbridgemode();
	
////===============检查开关动作是否正常，若出现异常则进行修正===============////
	//if (statechange[0]!=0)
	//{
	//	checkbridgemode();
	//	if (bridgeMode_1==0)
	//	{
	//		bridgeMode=bridgeMode_tmp+1;
	//		}
	//	if (bridgeMode_1!=0)
	//	{
	//		bridgeMode=bridgeMode_1+1;
	//	}
	//	while(bridgeMode>6)
	//	{
	//		bridgeMode=bridgeMode-6;
	//	}
	//}

	//if (statechange[1]!=0)
	//{
	//	checkbridgemode();
	//	if (bridgeMode_1!=0&&bridgeMode==0)
	//	{
	//		bridgeMode_tmp=bridgeMode_1;
	//	}
	//}
////========================================================================////

	if(statechange[0]+statechange[1] == 0)
	{
		return 0;
	}
	else
	{
		//cout<<"time is:"<<time<<endl;	
		//for (int i=0;i<6;i++)
		//	{
		//		cout<<diodeState[i]<<'\t';
		//	}
		//	cout<<endl;
		//switchtime_1=switchtime;
		//switchtime=time;
		//if (switchtime-switchtime_1<6.0e-5&&switchtime>1.0e-4)
		//{
		//	bridgeMode=bridgeMode_1;
		//	modifydiodemode();
		//	return 0;
		//}
		//checkbridgemode();
		//cout<<"time is:"<<time<<endl;	
		//for (int i=0;i<6;i++)
		//{
		//	cout<<diodeState[i]<<'\t';
		//}
		//cout<<endl;
		//modifybridgemode();
		calculateNortonEquivalentResistance(time);
		return 1;
	}

	////===============================================================////
	//int statechange[2] = {0,0};
	////double tmp_ratio;
	//for (int i=0;i<6;i++)
	//{
	//	if(diodeState[i]==0 && diodeVoltage[i]>forwardVoltageDrop)
	//	{
	//		diodeState[i] = 1;
	//		statechange[0] += 1;
	//		//tmp_ratio=diodeCurrent_1[i]/(diodeCurrent_1[i]-diodeCurrent[i]);
	//		//if (tmp_ratio<ratio_1)
	//		//{
	//		//	ratio_1=tmp_ratio;
	//		//}
	//	}
	//	else if(diodeState[i]==1 && diodeCurrent[i]<0)
	//	{
	//		diodeState[i] = 0;
	//		statechange[1] += 1;
	//		//tmp_ratio=diodeCurrent_1[i]/(diodeCurrent_1[i]-diodeCurrent[i]);
	//		//if (tmp_ratio<ratio_1)
	//		//{
	//		//	ratio_1=tmp_ratio;
	//		//}
	//		//ratio_1 = diodeCurrent_1[i]/(diodeCurrent_1[i]-diodeCurrent[i]);
	//	}
	//}
	//calculateNortonEquivalentResistance();
	//if(statechange[1] == 0)
	//{
	//	on_off = 1;
	//}
	//else
	//	on_off = 0;
	//if(statechange[0]+statechange[1] == 0)
	//{
	//	return 0;
	//}
	//else
	//	//return 1;
	//{
	//	cout<<time<<endl;
	//	for (int i=0;i<6;i++)
	//	{
	//		cout<<diodeState[i]<<'\t';
	//	}
	//	cout<<endl;
	//	return 1;
	//}
}

void ThreePhaseRectifierBridge::modifybridgemode()
{
	if (statechange[0]!=0)
	{
		//checkbridgemode();
		if (bridgeMode_1==0)
		{
			bridgeMode=bridgeMode_tmp+1;
		}
		if (bridgeMode_1!=0)
		{
			bridgeMode=bridgeMode_1+1;
		}
		while(bridgeMode>6)
		{
			bridgeMode=bridgeMode-6;
		}
	}
	if (bridgeMode_1!=0&&bridgeMode==0)
	{
		bridgeMode_tmp=bridgeMode_1;
	}
	modifydiodemode();
}
void ThreePhaseRectifierBridge::modifydiodemode()
{	
	for (int i=0;i<6;i++)
	{
		diodeState[i]=0;
	}
	if(bridgeMode==1)
	{
		diodeState[0]=1;
		diodeState[1]=1;
	}
	if(bridgeMode==2)
	{
		diodeState[1]=1;
		diodeState[2]=1;
	}
	if(bridgeMode==3)
	{
		diodeState[2]=1;
		diodeState[3]=1;
	}
	if(bridgeMode==4)
	{
		diodeState[3]=1;
		diodeState[4]=1;
	}
	if(bridgeMode==5)
	{
		diodeState[4]=1;
		diodeState[5]=1;
	}
	if(bridgeMode==6)
	{
		diodeState[5]=1;
		diodeState[0]=1;
	}
}
	
void ThreePhaseRectifierBridge::checkbridgemode()
{
	int sumdiodestate=0;
	for (int i=0;i<6;i++)
	{
		sumdiodestate=sumdiodestate+diodeState[i];
	}

	if (sumdiodestate==0)
	{
		bridgeMode=0;
	}
	else if (sumdiodestate==2)
	{
		if (diodeState[0]==1&&diodeState[1]==1)
		{
			bridgeMode=1;
		}
		if (diodeState[1]==1&&diodeState[2]==1)
		{
			bridgeMode=2;
		}
		if (diodeState[2]==1&&diodeState[3]==1)
		{
			bridgeMode=3;
		}
		if (diodeState[3]==1&&diodeState[4]==1)
		{
			bridgeMode=4;
		}
		if (diodeState[4]==1&&diodeState[5]==1)
		{
			bridgeMode=5;
		}
		if (diodeState[5]==1&&diodeState[0]==1)
		{
			bridgeMode=6;
		}
	}
	else
	{
		modifybridgemode();
		//	cout<<"123"<<endl;
		//	if (statechange[0]!=0)
		//	{
		//		//checkbridgemode();
		//		if (bridgeMode_1==0)
		//		{
		//			bridgeMode=bridgeMode_tmp+1;
		//		}
		//		if (bridgeMode_1!=0)
		//		{
		//			bridgeMode=bridgeMode_1+1;
		//		}
		//		while(bridgeMode>6)
		//		{
		//			bridgeMode=bridgeMode-6;
		//		}
	}

	//	if (statechange[1]!=0)
	//	{
	//		//checkbridgemode();
	//		if (bridgeMode_1!=0&&bridgeMode==0)
	//		{
	//			bridgeMode_tmp=bridgeMode_1;
	//		}
	//	}
	//}
	//modifybridgemode();
}

void ThreePhaseRectifierBridge::interpolate(double ratio)
{//给定比值，对支路的电压电流进行插值
	//for (int i=0;i<6;i++)
	//{
	//	diodeCurrent[i] = (1-ratio)*diodeCurrent_1[i] + ratio*diodeCurrent[i];
	//	diodeVoltage[i] = (1-ratio)*diodeVoltage_1[i] + ratio*diodeVoltage[i];
	//}
}

//int ThreePhaseRectifierBridge::getState()  //????
//{
//	/*return on_off;*/
//}

void ThreePhaseRectifierBridge::calculatebranchCurrentArray()
{
	branchCurrentArray[0] = diodeCurrent[0] - diodeCurrent[3];
	branchCurrentArray[1] = diodeCurrent[2] - diodeCurrent[5];
	branchCurrentArray[2] = diodeCurrent[4] - diodeCurrent[1];
	branchCurrentArray[3] = diodeCurrent[0] + diodeCurrent[2] + diodeCurrent[4];
}

//double ThreePhaseRectifierBridge::getratio()
//{
//	return ratio_1;
//}

//void ThreePhaseRectifierBridge:: printInfo()
//{
//
//
//	cout<<"nortonEquivalentResistance info:"<<endl;
//	for (int i=0;i<6;i++)
//	{
//		cout<<nortonEquivalentResistance[i]<<endl;
//	}
//
//}