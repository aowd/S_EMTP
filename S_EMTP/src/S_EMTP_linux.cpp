#ifdef linux

#include "EMTP.h"

#include <sys/time.h>
#include <unistd.h>


void S_EMTP_linux()
{
	struct timeval tpstart,tpend; 
	float timeuse; 

	EMTP emtp;
	emtp.specifySystem();//ϵͳ���뺯����ʹ��Ĭ�ϲ���������ʱ��60ms�����沽��50us��

//������ѭ��
		emtp.counter=1;//������
		emtp.curTime=0.0;//��ǰʱ��

		emtp.initializeSystem(0);//��ʼ��ϵͳ��0--0������1--PSCAD��ֵ����
		emtp.lu_conductanceMatrix= new TDecompLU(emtp.conductanceMatrix);
		emtp.formConductanceMatrix();

		//Update current time and counter
			emtp.curTime= emtp.curTime +emtp.deltaT;
			emtp.counter =emtp.counter + 1;
			emtp.timeArray(emtp.counter) = emtp.curTime;
			emtp.checkSwitch();// Check switch state			
			emtp.calculateBranchNortonEquivalentCurrent();//����֧·ŵ�ٵ�Ч����
			emtp.formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������
//===================��ʱ��ʼ===================
		int nstep=0;		
		gettimeofday(&tpstart,NULL);
		
		while(emtp.curTime<emtp.finishTime-emtp.deltaT/2)
		{
			nstep++;
			emtp.solveNodeVoltageEquation();//���ڵ��ѹ����

			emtp.saveNodeVoltage();//����ڵ��ѹ����,�����������ȫ�ֵ�ѹ�����У����֧·�޹ء�

			emtp.calculateBranchVoltage();//����֧·��ѹ������

			emtp.calculateBranchCurrent();//����֧·����

			emtp.saveBranchCurrent();//����֧·����
		
			emtp.curTime= emtp.curTime +emtp.deltaT;
			emtp.counter =emtp.counter + 1;
			emtp.timeArray(emtp.counter) = emtp.curTime;
			emtp.checkSwitch();// Check switch state
			emtp.calculateBranchNortonEquivalentCurrent();//����֧·ŵ�ٵ�Ч����
			emtp.formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������

		}

		emtp.solveNodeVoltageEquation();//���ڵ��ѹ����

		emtp.saveNodeVoltage();//����ڵ��ѹ����,�����������ȫ�ֵ�ѹ�����У����֧·�޹ء�

		emtp.calculateBranchVoltage();//����֧·��ѹ������

		emtp.calculateBranchCurrent();//����֧·����

		emtp.saveBranchCurrent();//����֧·����
		emtp.checkSwitch();// Check switch state

		gettimeofday(&tpend,NULL); 
//===================��ʱ����===================		
		emtp.saveResult();//ϵͳ�������������������2ά����ĸ�ʽ�洢���ı��ļ�result.dat��
//		emtp.checkResult();
         

	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+ tpend.tv_usec-tpstart.tv_usec; 
//	timeuse/=1000;
	std::cout<<"case "<<emtp.caseNumber<<" average step time cost "<<timeuse/nstep<<"us!"<<endl;

	delete emtp.lu_conductanceMatrix;//ÿһ��case����new������lu����Ҫɾ����
}

#endif