#ifdef linux

#include "EMTP.h"

#include <sys/time.h>
#include <unistd.h>


void S_EMTP_linux()
{
	struct timeval tpstart,tpend; 
	float timeuse; 

	EMTP emtp;
	emtp.specifySystem();//系统输入函数：使用默认参数，仿真时间60ms，仿真步长50us。

//进入主循环
		emtp.counter=1;//计数器
		emtp.curTime=0.0;//当前时间

		emtp.initializeSystem(0);//初始化系统，0--0启动；1--PSCAD初值启动
		emtp.lu_conductanceMatrix= new TDecompLU(emtp.conductanceMatrix);
		emtp.formConductanceMatrix();

		//Update current time and counter
			emtp.curTime= emtp.curTime +emtp.deltaT;
			emtp.counter =emtp.counter + 1;
			emtp.timeArray(emtp.counter) = emtp.curTime;
			emtp.checkSwitch();// Check switch state			
			emtp.calculateBranchNortonEquivalentCurrent();//计算支路诺顿等效电流
			emtp.formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量
//===================计时开始===================
		int nstep=0;		
		gettimeofday(&tpstart,NULL);
		
		while(emtp.curTime<emtp.finishTime-emtp.deltaT/2)
		{
			nstep++;
			emtp.solveNodeVoltageEquation();//求解节点电压方程

			emtp.saveNodeVoltage();//保存节点电压向量,将求解结果存在全局电压矩阵中，与各支路无关。

			emtp.calculateBranchVoltage();//计算支路电压并保存

			emtp.calculateBranchCurrent();//计算支路电流

			emtp.saveBranchCurrent();//保存支路电流
		
			emtp.curTime= emtp.curTime +emtp.deltaT;
			emtp.counter =emtp.counter + 1;
			emtp.timeArray(emtp.counter) = emtp.curTime;
			emtp.checkSwitch();// Check switch state
			emtp.calculateBranchNortonEquivalentCurrent();//计算支路诺顿等效电流
			emtp.formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量

		}

		emtp.solveNodeVoltageEquation();//求解节点电压方程

		emtp.saveNodeVoltage();//保存节点电压向量,将求解结果存在全局电压矩阵中，与各支路无关。

		emtp.calculateBranchVoltage();//计算支路电压并保存

		emtp.calculateBranchCurrent();//计算支路电流

		emtp.saveBranchCurrent();//保存支路电流
		emtp.checkSwitch();// Check switch state

		gettimeofday(&tpend,NULL); 
//===================计时结束===================		
		emtp.saveResult();//系统输出函数：将计算结果按2维数组的格式存储在文本文件result.dat中
//		emtp.checkResult();
         

	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+ tpend.tv_usec-tpstart.tv_usec; 
//	timeuse/=1000;
	std::cout<<"case "<<emtp.caseNumber<<" average step time cost "<<timeuse/nstep<<"us!"<<endl;

	delete emtp.lu_conductanceMatrix;//每一个case里面new出来的lu对象都要删除掉
}

#endif