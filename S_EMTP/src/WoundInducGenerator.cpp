#include "WoundInducGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                ���캯������������                                    */
/************************************************************************/
WoundInducGenerator::WoundInducGenerator(int id, int firstNode,int lastNode, int control,double Vw)
{
	this->id = id;
	isUserDef = 0;
	need_NEC=1;
	//���캯����a��b��cΪ�ڵ���
	nPort=6;
	nodeNumber=new int[nPort];
	type=18;

	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	Wbase = 2*PI*50;//����ٶȣ�����ֵ,f=60Hz
	w=1;
	
	//dq0����ϵ�µĵ������,����ֵ
	Rs = 0.048;
	Rr = 0.018;
	//dq0����ϵ�µĵ翹����,����ֵ
	Xm = 3.8;
	Xls = 0.075;
	Xlr = 0.12;
	
	//���Ӳ��ѹ���������迹�Ļ�ֵ
	Vbase=690;
	Sbase=2e6;
	Vbase = Vbase/sqrt(3.0);//V
	Ibase = Sbase/3/Vbase;//A
	Zbase = Vbase/Ibase;
	
	//��ת������������
	turnratio=0.3;

	this->control = control;

	//��е�����
	H = 3;
	D = 0;
	TL =-1.0;//����ת��

	//�������
	p_air=1.225;
	rArea=3177;
	GR=90.5;
	beta=0;
	GE=0.979;
	vWind=Vw-0.1*(id-1);

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/

	//�ɿ˱任�����������
	Ps.ResizeTo(1,2,1,3);
	invPs.ResizeTo(1,3,1,2);
	Pr.ResizeTo(1,2,1,3);
	invPr.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,4,1,4);
	BB.ResizeTo(1,4,1,4);
	AA_inv.ResizeTo(1,4,1,4);

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ�����������0����ʼ��ʱ������ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vsabc[i-1] = 0;
		Vsabc_1[i-1] = 0;
		Vsabc_2[i-1] = 0;
		Isabc[i-1] = 0;
		Isabc_1[i-1] = 0;
		Isabc_2[i-1] = 0;
		Isrc[i-1]=0;
		Isrc_1[i-1]=0;
		Isrc_2[i-1]=0;
	}
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=2;i++)
	{
		Vsdq[i-1] = 0;
		Isdq[i-1] = 0;
	}

	//abc����ϵ��ת����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vrabc[i-1] = 0;
		Vrabc_1[i-1] = 0;
		Vrabc_2[i-1] = 0;
		Irabc[i-1] = 0;
		Irabc_1[i-1] = 0;
		Irabc_2[i-1] = 0;
		Irrc[i-1]=0;
		Irrc_1[i-1]=0;
		Irrc_2[i-1]=0;
	}
	//dq0����ϵ��ת����Ȧ�ĵ�ѹ����
	for (int i=1;i<=2;i++)
	{
		Vrdq[i-1] = 0;
		Irdq[i-1] = 0;
	}

	/*  ��е��������  */
	if (control==1)
	{
		wr = 1.16715;//ת��ת��
		wr_his = 1.16715;//wr����ʷֵ
	} 
	else
	{
		wr = 0;//ת��ת��
		wr_his = 0;//wr����ʷֵ
	}
	curAngle_s= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	curAngle_r= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Te=0;
	Te_his=0;


	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	for (int i=0;i<2;i++)
	{
		nortonEquivalentConductance[i] = 0;
		nortonEquivalentConductance_pu[i]=0;
		RCnortonEquivalentConductance[i]=0;
	}

	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
		RCnortonEquivalentCurrent[i-1]=0;
		RCnortonEquivalentCurrent_1[i-1]=0;
		RCnortonEquivalentCurrent_2[i-1]=0;
	}
	for(int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq[i-1] = 0;//dq0�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_dq_his[i-1] = 0;
	}

	/************************************************************************/
	/*              ��EMTP�ӿ�����ı���                                    */
	/************************************************************************/
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}
#ifdef WIND_VELOCITY_DATA_INPUT
	WindVelocityCounter=0;
	WindVelocityCounter_bak=0;
#endif
}
WoundInducGenerator::~WoundInducGenerator()
{//��������

}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/

void WoundInducGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time){ptr+=6;}

void WoundInducGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
{//�ӽڵ��ѹ�����ж�֧·�ڵ��ѹ

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
		Vsabc_his2[i] = Vsabc_his[i];
		Vsabc_his[i] = Vsabc[i];
	}

	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i]==0)
		{
			Vsabc[i]=0;
		}
		else
		{
			Vsabc[i]= nodeVoltageArray(nodeNumber[i]);
		}	
	}

	for (int i=0;i<3;i++)
	{
		Vrabc_1[i] =Vrabc[i];
		Vrabc_his2[i] = Vrabc_his[i];
		Vrabc_his[i] = Vrabc[i];
	}


	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i+3]==0)
		{
			Vrabc[i]=0;
		}
		else
		{
			Vrabc[i]= nodeVoltageArray(nodeNumber[i+3]);
		}	
	}
}
void WoundInducGenerator::calculateBranchVoltage()
{//����֧·��ѹ
}

void WoundInducGenerator::calculateBranchCurrent()
{//����֧·����
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
#ifdef S_MACHINE_WITH_COMPENSATE_RESISTANCE
		Isabc[i] = Vsabc[i]*nortonEquivalentConductance[0]+nortonEquivalentCurrent[i];
#endif
#ifdef S_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
		Isabc[i] = nortonEquivalentCurrent[i];
#endif
		Irabc_1[i] = Irabc[i];
#ifdef R_MACHINE_WITH_COMPENSATE_RESISTANCE
		Irabc[i] = Vrabc[i]*nortonEquivalentConductance[1]+nortonEquivalentCurrent[i+3];
#endif
#ifdef R_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
		Irabc[i] = nortonEquivalentCurrent[i+3];
#endif
		
		Isrc_1[i] = Isrc[i];
		Isrc[i] = Vsabc[i]*RCnortonEquivalentConductance[0]+RCnortonEquivalentCurrent[i];
		Irrc_1[i] = Irrc[i];
		Irrc[i] = Vrabc[i]*RCnortonEquivalentConductance[1]+RCnortonEquivalentCurrent[i+3];
	}
}

void WoundInducGenerator::calculateNortonEquivalentCurrent(double time)
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	
	calculateDqResults();//����dq0����ϵ�µĽ��
	calculateOmega();
	calculateCoefficientMatrix();

	////����ABC����ϵ�µĵ�ѹ��ֵԤ��
	double Vsdq_pre[2],Vrdq_pre[2];
	for (int i=0;i<2;i++)
	{
		Vsdq_pre[i] = Ps(i+1,1)*(1.25*Vsabc[0]+0.5*Vsabc_his[0]-0.75*Vsabc_his2[0])+Ps(i+1,2)*(1.25*Vsabc[1]+0.5*Vsabc_his[1]-0.75*Vsabc_his2[1])
			+Ps(i+1,3)*(1.25*Vsabc[2]+0.5*Vsabc_his[2]-0.75*Vsabc_his2[2]);
		Vrdq_pre[i] = Pr(i+1,1)*(1.25*Vrabc[0]+0.5*Vrabc_his[0]-0.75*Vrabc_his2[0])+Pr(i+1,2)*(1.25*Vrabc[1]+0.5*Vrabc_his[1]-0.75*Vrabc_his2[1])
			+Pr(i+1,3)*(1.25*Vrabc[2]+0.5*Vrabc_his[2]-0.75*Vrabc_his2[2]);
		Vsdq_pre[i] = Vsdq_pre[i]/Vbase;
		Vrdq_pre[i] = Vrdq_pre[i]/Vbase*turnratio;
	}

	TMatrixD CC=AA_inv*BB;
	//����dq0�����µ���ʷ������
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq_his[i-1] = nortonEquivalentCurrent_dq[i-1];
		
		//������֧·����
#ifdef S_MACHINE_WITHOUT_INT
		//�����õ�ѹ��ֵԤ��
		nortonEquivalentCurrent_dq[i-1] =-(CC(i,1)*Isdq[0]+CC(i,2)*Isdq[1]+CC(i,3)*Irdq[0]+CC(i,4)*Irdq[1])
		+ 2*(AA_inv(i,1)*Vsdq[0]+AA_inv(i,2)*Vsdq[1]+AA_inv(i,3)*Vrdq[0]+AA_inv(i,4)*Vrdq[1]);
#endif

#ifdef S_ABC_MACHINE_WITH_INT
		//����ABC���ѹ��ֵԤ��
		nortonEquivalentCurrent_dq[i-1] =-(CC(i,1)*Isdq[0]+CC(i,2)*Isdq[1]+CC(i,3)*Irdq[0]+CC(i,4)*Irdq[1])
			+ AA_inv(i,1)*Vsdq_pre[0]+AA_inv(i,2)*Vsdq_pre[1]+AA_inv(i,3)*Vrdq_pre[0]+AA_inv(i,4)*Vrdq_pre[1]
			+AA_inv(i,1)*Vsdq[0]+AA_inv(i,2)*Vsdq[1]+AA_inv(i,3)*Vrdq[0]+AA_inv(i,4)*Vrdq[1];
#endif

#ifdef S_DQ_MACHINE_WITH_INT
		////����DQ���ѹ��ֵԤ��
		nortonEquivalentCurrent_dq[i-1] =-(CC(i,1)*Isdq[0]+CC(i,2)*Isdq[1]+CC(i,3)*Irdq[0]+CC(i,4)*Irdq[1])
			+ AA_inv(i,1)*(1.25*Vsdq[0]+0.5*Vsdq_his[0]-0.75*Vsdq_his2[0])+AA_inv(i,2)*(1.25*Vsdq[1]+0.5*Vsdq_his[1]-0.75*Vsdq_his2[1])
			+AA_inv(i,3)*(1.25*Vrdq[0]+0.5*Vrdq_his[0]-0.75*Vrdq_his2[0])+AA_inv(i,4)*(1.25*Vrdq[1]+0.5*Vrdq_his[1]-0.75*Vrdq_his2[1])
			+AA_inv(i,1)*Vsdq[0]+AA_inv(i,2)*Vsdq[1]+AA_inv(i,3)*Vrdq[0]+AA_inv(i,4)*Vrdq[1];
#endif

		//����DQ�Ჹ������
#ifdef S_DQ_COMPENSATE_WITHOUT_INT
		nortonEquivalentCurrent_dq[i-1] =nortonEquivalentCurrent_dq[i-1]-Vsdq[i-1]*nortonEquivalentConductance_pu[0];
#endif
#ifdef S_DQ_COMPENSATE_WITH_INT
		nortonEquivalentCurrent_dq[i-1] =nortonEquivalentCurrent_dq[i-1]-(1.25*Vsdq[i-1]+0.5*Vsdq_his[i-1]-0.75*Vsdq_his2[i-1])*nortonEquivalentConductance_pu[0];
		//nortonEquivalentCurrent_dq[i-1] =nortonEquivalentCurrent_dq[i-1]-(4.0/3.0*Vsdq[i-1]+1.0/3.0*Vsdq_his[i-1]-2.0/3.0*Vsdq_his2[i-1])*nortonEquivalentConductance_pu[0];
#endif

		nortonEquivalentCurrent_dq_his[i+1] = nortonEquivalentCurrent_dq[i+1];
		
		//������֧·����
#ifdef R_MACHINE_WITHOUT_INT	//�����õ�ѹ��ֵԤ��
		nortonEquivalentCurrent_dq[i+1] =-(CC(i+2,1)*Isdq[0]+CC(i+2,2)*Isdq[1]+CC(i+2,3)*Irdq[0]+CC(i+2,4)*Irdq[1])
			+ 2*(AA_inv(i+2,1)*Vsdq[0]+AA_inv(i+2,2)*Vsdq[1]+AA_inv(i+2,3)*Vrdq[0]+AA_inv(i+2,4)*Vrdq[1]);
#endif
#ifdef R_ABC_MACHINE_WITH_INT		//����ABC���ѹ��ֵԤ��
		nortonEquivalentCurrent_dq[i+1] =-(CC(i+2,1)*Isdq[0]+CC(i+2,2)*Isdq[1]+CC(i+2,3)*Irdq[0]+CC(i+2,4)*Irdq[1])
			+ AA_inv(i+2,1)*Vsdq_pre[0]+AA_inv(i+2,2)*Vsdq_pre[1]+AA_inv(i+2,3)*Vrdq_pre[0]+AA_inv(i+2,4)*Vrdq_pre[1]
			+AA_inv(i+2,1)*Vsdq[0]+AA_inv(i+2,2)*Vsdq[1]+AA_inv(i+2,3)*Vrdq[0]+AA_inv(i+2,4)*Vrdq[1];
#endif
#ifdef R_DQ_MACHINE_WITH_INT 	 //����DQ���ѹ��ֵ
		nortonEquivalentCurrent_dq[i+1] =-(CC(i+2,1)*Isdq[0]+CC(i+2,2)*Isdq[1]+CC(i+2,3)*Irdq[0]+CC(i+2,4)*Irdq[1])
			+ AA_inv(i+2,1)*(1.25*Vsdq[0]+0.5*Vsdq_his[0]-0.75*Vsdq_his2[0])+AA_inv(i+2,2)*(1.25*Vsdq[1]+0.5*Vsdq_his[1]-0.75*Vsdq_his2[1])
			+AA_inv(i+2,3)*(1.25*Vrdq[0]+0.5*Vrdq_his[0]-0.75*Vrdq_his2[0])+AA_inv(i+2,4)*(1.25*Vrdq[1]+0.5*Vrdq_his[1]-0.75*Vrdq_his2[1])
			+AA_inv(i+2,1)*Vsdq[0]+AA_inv(i+2,2)*Vsdq[1]+AA_inv(i+2,3)*Vrdq[0]+AA_inv(i+2,4)*Vrdq[1];
#endif

		//����DQ�Ჹ������
#ifdef R_DQ_COMPENSATE_WITHOUT_INT
		nortonEquivalentCurrent_dq[i+1] =nortonEquivalentCurrent_dq[i+1]-Vrdq[i-1]*nortonEquivalentConductance_pu[1];
#endif
#ifdef R_DQ_COMPENSATE_WITH_INT
		nortonEquivalentCurrent_dq[i+1] =nortonEquivalentCurrent_dq[i+1]-(1.25*Vrdq[i-1]+0.5*Vrdq_his[i-1]-0.75*Vrdq_his2[i-1])*nortonEquivalentConductance_pu[1];
		//nortonEquivalentCurrent_dq[i+1] =nortonEquivalentCurrent_dq[i+1]-(4.0/3.0*Vrdq[i-1]+1.0/3.0*Vrdq_his[i-1]-2.0/3.0*Vrdq_his2[i-1])*nortonEquivalentConductance_pu[1];
#endif
	}

	//����park��任����
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//����abc�����µ���ʷ������
	for (int i=1;i<=6;i++)
	{
		if (i<=3)
		{
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq[0]+ invPs(i,2)*nortonEquivalentCurrent_dq[1];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;

#ifdef S_ABC_COMPENSATE_WITHOUT_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-Vsabc[i-1]*nortonEquivalentConductance[0];//ABC�Ჹ������
#endif
#ifdef S_ABC_COMPENSATE_WITH_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(1.25*Vsabc[i-1]+0.5*Vsabc_his[i-1]-0.75*Vsabc_his2[i-1])*nortonEquivalentConductance[0];//ABC�Ჹ������������ֵ
			//nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(4.0/3.0*Vsabc[i-1]+1.0/3.0*Vsabc_his[i-1]-2.0/3.0*Vsabc_his2[i-1])*nortonEquivalentConductance[0];//ABC�Ჹ������������ֵ
#endif
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq[2]+ invPr(i-3,2)*nortonEquivalentCurrent_dq[3];			
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;

#ifdef R_ABC_COMPENSATE_WITHOUT_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-Vrabc[i-4]*nortonEquivalentConductance[1];//ABC�Ჹ������
#endif
#ifdef R_ABC_COMPENSATE_WITH_INT
			nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(1.25*Vrabc[i-4]+0.5*Vrabc_his[i-4]-0.75*Vrabc_his2[i-4])*nortonEquivalentConductance[1];//ABC�Ჹ������������ֵ
			//nortonEquivalentCurrent[i-1]= nortonEquivalentCurrent[i-1]-(4.0/3.0*Vrabc[i-4]+1.0/3.0*Vrabc_his[i-4]-2.0/3.0*Vrabc_his2[i-4])*nortonEquivalentConductance[1];//ABC�Ჹ������������ֵ
#endif
		}		
	}

	//����RC֧·��ŵ�ٵ�ֵ����
	double Pc[2];
	Pc[0]=deltaT/2/Crc[0];
	Pc[1]=deltaT/2/Crc[1];
	for (int i=0;i<3;i++)
	{
		RCnortonEquivalentCurrent[i] = (Rrc[0]-Pc[0])/(Rrc[0]+Pc[0])*Isrc[i] - 1/(Rrc[0]+Pc[0])*Vsabc[i];
		RCnortonEquivalentCurrent[i+3] = (Rrc[1]-Pc[1])/(Rrc[1]+Pc[1])*Irrc[i] - 1/(Rrc[1]+Pc[1])*Vrabc[i];
	}
}

void WoundInducGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//�γɽڵ�ŵ�ٵ�Ч��������
	int N;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) = nodeNortonEquivalentCurrentArray(N)- nortonEquivalentCurrent[i] - RCnortonEquivalentCurrent[i];
	}
}
void WoundInducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//�γɽڵ㵼����
	int N;
	double tmpd;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		if (i<3)
		{
#ifdef S_MACHINE_WITH_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[0]+RCnortonEquivalentConductance[0];
#endif
#ifdef S_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+RCnortonEquivalentConductance[0];
#endif
		} 
		else
		{
#ifdef R_MACHINE_WITH_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[1]+RCnortonEquivalentConductance[1];
#endif
#ifdef R_MACHINE_WITHOUT_COMPENSATE_RESISTANCE
			conductanceMatrix(N,N) = tmpd+RCnortonEquivalentConductance[1];
#endif
		}
	}
}

void WoundInducGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//����֧·����
	branchCurrentMatrix(counter,ptr) = Isabc[0]+Isrc[0];
	branchCurrentMatrix(counter,ptr+1) = Isabc[1]+Isrc[1];
	branchCurrentMatrix(counter,ptr+2) = Isabc[2]+Isrc[2];
	ptr+=3;
	branchCurrentMatrix(counter,ptr) = Irabc[0]+Irrc[0];
	branchCurrentMatrix(counter,ptr+1) = Irabc[1]+Irrc[1];
	branchCurrentMatrix(counter,ptr+2) = Irabc[2]+Irrc[2];
	ptr+=3;
}

void WoundInducGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//����֧·����
	branchCurrentMatrix_1[counter-1][ptr] = Isabc[0]+Isrc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Isabc[1]+Isrc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Isabc[2]+Isrc[2];
	ptr+=3;
	branchCurrentMatrix_1[counter-1][ptr] = Irabc[0]+Irrc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Irabc[1]+Irrc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Irabc[2]+Irrc[2];
	ptr+=3;
}
void WoundInducGenerator::interpolate(double ratio)
{//������ֵ����֧·�ĵ�ѹ�������в�ֵ

	//��ѹ����������ֵ
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Isrc[k] = (1-ratio)*Isrc_1[k] + ratio*Isrc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
		Irrc[k] = (1-ratio)*Irrc_1[k] + ratio*Irrc[k];
	}

	for (int k=0;k<2;k++)
	{
		Vsdq[k] = (1-ratio)*Vsdq_his[k] + ratio*Vsdq[k];
		Isdq[k] = (1-ratio)*Isdq_his[k] + ratio*Isdq[k];
		Vrdq[k] = (1-ratio)*Vrdq_his[k] + ratio*Vrdq[k];
		Irdq[k] = (1-ratio)*Irdq_his[k] + ratio*Irdq[k];
	}

	for (int k=0;k<4;k++)
	{
		nortonEquivalentCurrent_dq[k] = (1-ratio)*nortonEquivalentCurrent_dq_his[k] + ratio*nortonEquivalentCurrent_dq[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent[k] = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		RCnortonEquivalentCurrent[k] = (1-ratio)*RCnortonEquivalentCurrent_1[k] + ratio*RCnortonEquivalentCurrent[k];
	}

	curAngle_r = (1-ratio)*curAngle_r_1 + ratio*curAngle_r;
	curAngle_s = (1-ratio)*curAngle_s_1 + ratio*curAngle_s;
	wr = (1-ratio)*wr_his + ratio*wr;

 	//ŵ�ٵ�ֵ������ֵ
	//double temp1,temp2;
	//for (int k=0;k<6;k++)
	//{
	//	temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	//	temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

	//	nortonEquivalentCurrent_1[k] = temp1;
	//	nortonEquivalentCurrent[k] = temp2;
	//}
}

void WoundInducGenerator::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for (int i=0;i<3;i++)
		{
			Vsabc_2[i] = Vsabc_1[i];
			Isabc_2[i] = Isabc_1[i];
			Isrc_2[i] = Isrc_1[i];
			Vrabc_2[i] = Vrabc_1[i];
			Irabc_2[i] = Irabc_1[i];
			Irrc_2[i] = Irrc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_2[i+3] = nortonEquivalentCurrent_1[i+3];
			RCnortonEquivalentCurrent_2[i] = RCnortonEquivalentCurrent_1[i];
			RCnortonEquivalentCurrent_2[i+3] = RCnortonEquivalentCurrent_1[i+3];
		}
		break;
	case 2://��_2��������ֵ����_1������
		for (int i=0;i<3;i++)
		{
			Vsabc_1[i] = Vsabc_2[i];
			Isabc_1[i] = Isabc_2[i];
			Isrc_1[i] = Isrc_2[i];
			Vrabc_1[i] = Vrabc_2[i];
			Irabc_1[i] = Irabc_2[i];
			Irrc_1[i] = Irrc_2[i];
		}
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
			RCnortonEquivalentCurrent[i] = RCnortonEquivalentCurrent_1[i];
			RCnortonEquivalentCurrent_1[i] = RCnortonEquivalentCurrent_2[i];
		}
		break;
	case 3://ŵ�ٵ�ֵ�����洢����_1��������ֵ����_2������
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			RCnortonEquivalentCurrent[i] = RCnortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"��������ȷ�ĸ������ͱ�ţ�"<<endl;
		exit(1);
	}
}

/************************************************************************/
/*                ������ڲ����õĺ���                                  */
/************************************************************************/
//����ϵ������AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs4,Gr1...Gr3
//ͬʱ����nortonEquivalentConductance
void WoundInducGenerator::calculateCoefficientMatrix()
{

	//����AA���󣬼��������
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w*(Xm+Xls);
	AA(1,3) = T*Xm;
	AA(1,4) = -w*Xm;

	AA(2,1) = w*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,3) = w*Xm;
	AA(2,4) = T*Xm;

	AA(3,1) = T*Xm;
	AA(3,2) = (wr-w)*Xm;
	AA(3,3) = Rr+T*(Xm+Xlr);
	AA(3,4) = (wr-w)*(Xm+Xlr);

	AA(4,1) = (w-wr)*Xm;
	AA(4,2) = T*Xm;
	AA(4,3) = (w-wr)*(Xm+Xlr);
	AA(4,4) = Rr+T*(Xm+Xlr);

	//����BB����
	BB(1,1) = Rs-T*(Xm+Xls);
	BB(1,2) = -w*(Xm+Xls);
	BB(1,3) = -T*Xm;
	BB(1,4) = -w*Xm;

	BB(2,1) = w*(Xm+Xls);
	BB(2,2) = Rs-T*(Xm+Xls);
	BB(2,3) = w*Xm;
	BB(2,4) = -T*Xm;

	BB(3,1) = -T*Xm;
	BB(3,2) = (wr-w)*Xm;
	BB(3,3) = Rr-T*(Xm+Xlr);
	BB(3,4) = (wr-w)*(Xm+Xlr);

	BB(4,1) = (w-wr)*Xm;
	BB(4,2) = -T*Xm;
	BB(4,3) = (w-wr)*(Xm+Xlr);
	BB(4,4) = Rr-T*(Xm+Xlr);

	AA_inv = AA;
	AA_inv.Invert();
}

//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
void WoundInducGenerator::calculateDqResults()
{
	double Vsabc_pu[3];
	double Isabc_pu[3];
	double Vrabc_pu[3];
	double Irabc_pu[3];

	//��abc�����µĶ��ӵ�ѹ����ת��Ϊ����ֵ
	for (int i=0;i<3;i++)
	{
		Vsabc_pu[i] = Vsabc[i]/Vbase;
		Isabc_pu[i] = Isabc[i]/Ibase;
		Vrabc_pu[i] = Vrabc[i]/Vbase*turnratio;
		Irabc_pu[i] = Irabc[i]/Ibase/turnratio;
	}

	//����park�任����
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	//�����ѹ
	for (int i=0;i<2;i++)
	{
		Vsdq_his2[i] = Vsdq_his[i];
		Vsdq_his[i] = Vsdq[i];
		Vsdq[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq_his2[i] = Vrdq_his[i];
		Vrdq_his[i] = Vrdq[i];
		Vrdq[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//�������
	for (int i=0;i<2;i++)
	{
		Isdq_his2[i] = Isdq_his[i];
		Isdq_his[i] = Isdq[i];
		Isdq[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq_his2[i] = Irdq_his[i];
		Irdq_his[i] = Irdq[i];
		Irdq[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}
}
//����е���̣�����ת��ת��
//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
void WoundInducGenerator::calculateOmega()
{
	//���÷������϶���ʽ
	if (control==2)
	{
		//���������ת��
		calculateTwind();
		//������ת��Ϊ�϶�ת��
		TL=-Twind;
	}
	wr_his = wr;
	curAngle_r_1= curAngle_r;
	curAngle_s_1=curAngle_s;

	if (control==0||control==2)
	{
		Te_his=Te;
		Te=0.5*Xm*(Irdq[0]*Isdq[1]-Irdq[1]*Isdq[0]);
		wr = ((Te+Te_his)/2-TL)*deltaT/2/H+wr_his;
		curAngle_r = curAngle_r + deltaT*Wbase*(w-(wr_his+wr)/2);
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(w-wr)*Wbase;
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
}

//�����ɿ˱任�뷴�任����
void WoundInducGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void WoundInducGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}


// ƽ��ģ����أ�xuyin��20121224
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
void WoundInducGenerator::storeInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc_bak[k] = Vsabc[k];
		Vsabc_his_bak[k] = Vsabc_his[k];
		Vsabc_his2_bak[k] = Vsabc_his2[k];
		Isabc_bak[k] = Isabc[k];
		Isrc_bak[k] = Isrc[k];
		Vrabc_bak[k] = Vrabc[k];		
		Vrabc_his_bak[k] = Vrabc_his[k];
		Vrabc_his2_bak[k] = Vrabc_his2[k];
		Irabc_bak[k] = Irabc[k];
		Irrc_bak[k] = Irrc[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq_bak[k] = Vsdq[k];
		Vsdq_his_bak[k] = Vsdq_his[k];
		Vsdq_his2_bak[k] = Vsdq_his2[k];
		Isdq_bak[k] = Isdq[k];
		Isdq_his_bak[k] = Isdq_his[k];
		Isdq_his2_bak[k] = Isdq_his2[k];
		Vrdq_bak[k] = Vrdq[k];
		Vrdq_his_bak[k] = Vrdq_his[k];
		Vrdq_his2_bak[k] = Vrdq_his2[k];
		Irdq_bak[k] = Irdq[k];
		Irdq_his_bak[k] = Irdq_his[k];
		Irdq_his2_bak[k] = Irdq_his2[k];
	}

	wr_bak = wr;
	wr_his_bak = wr_his;
	Te_bak = Te;
	Te_his_bak=Te_his;
	curAngle_s_bak = curAngle_s;
	curAngle_r_bak = curAngle_r;
	WindVelocityCounter_bak = WindVelocityCounter;
}

// ����У��ʱ�ָ��ڲ�����
void WoundInducGenerator::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc[k] = Vsabc_bak[k];
		Vsabc_his[k] = Vsabc_his_bak[k];
		Vsabc_his2[k] = Vsabc_his2_bak[k];
		Isabc[k] = Isabc_bak[k];
		Isrc[k] = Isrc_bak[k];
		Vrabc[k] = Vrabc_bak[k];
		Vrabc_his[k] = Vrabc_his_bak[k];
		Vrabc_his2[k] = Vrabc_his2_bak[k];
		Irabc[k] = Irabc_bak[k];
		Irrc[k] = Irrc_bak[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq[k] = Vsdq_bak[k];
		Vsdq_his[k] = Vsdq_his_bak[k];
		Vsdq_his2[k] = Vsdq_his2_bak[k];
		Isdq[k] = Isdq_bak[k];
		Isdq_his[k] = Isdq_his_bak[k];
		Isdq_his2[k] = Isdq_his2_bak[k];
		Vrdq[k] = Vrdq_bak[k];
		Vrdq_his[k] = Vrdq_his_bak[k];
		Vrdq_his2[k] = Vrdq_his2_bak[k];
		Irdq[k] = Irdq_bak[k];
		Irdq_his[k] = Irdq_his_bak[k];
		Irdq_his2[k] = Irdq_his2_bak[k];
	}

	wr = wr_bak;
	wr_his = wr_his_bak;
	Te = Te_bak;
	Te_his = Te_his_bak;
	curAngle_s = curAngle_s_bak;
	curAngle_r = curAngle_r_bak;
	WindVelocityCounter = WindVelocityCounter_bak;
}

void WoundInducGenerator::calculateNortonEquivalentResistance(double time)
{
	//���������翹
	double Lds,Ldr,nortonEquivalentResistance[2];
	Lds=(Xls+Xm*Xlr/(Xm+Xlr))*Zbase/Wbase;
	Ldr=Xlr*Zbase/Wbase;
	//���㲹������
	nortonEquivalentResistance[0]=2*Lds/deltaT;
	nortonEquivalentResistance[1]=2*Ldr/deltaT;
	nortonEquivalentConductance[0]=1/nortonEquivalentResistance[0];
	nortonEquivalentConductance[1]=1/nortonEquivalentResistance[1]*turnratio*turnratio;
	//����Ч������ֵ���㵽����ֵ
	nortonEquivalentConductance_pu[0] = nortonEquivalentConductance[0]*Zbase;
	nortonEquivalentConductance_pu[1] = nortonEquivalentConductance[1]*Zbase/turnratio/turnratio;
	
	T = 2.0/deltaT/Wbase;
	//����RC֧·����
	Rrc[0]=20*Zbase;
	Rrc[1]=Rrc[0]/turnratio/turnratio;
	Crc[0]=10*deltaT/Rrc[0];
	Crc[1]=10*deltaT/Rrc[1];
	//Crc[0]=10*(5e-6)/Rrc[0];
	//Crc[1]=10*(5e-6)/Rrc[1];
	RCnortonEquivalentConductance[0]=1/(Rrc[0]+deltaT/2/Crc[0]);
	RCnortonEquivalentConductance[1]=1/(Rrc[1]+deltaT/2/Crc[1]);
}

// �����ʼ��ר�ú���, xuyin, 20121226
void WoundInducGenerator::initializeGen(double** GenInitialMatrix, int& ptr)
{
	// ���˵�һ��Ϊʱ��t�⣬ÿ�������ӦGenInitialMatrix��15��
	// 1-3�У�Vsabc��4-6�У�Vrabc��7-9�У�Isabc��10-12�У�Irabc��13-15�У�Tm��Theta��wr
	// GenInitialMatrix�����У��ֱ��Ӧ��ǰʱ���Լ�ǰ����ʱ����ֵ������ʱ��˳��

	// ��������ֵ*1000
	for (int i=0; i<3; i++)
		for (int j=0; j<12; j++)
			GenInitialMatrix[i][ptr+j] *= 1000;

	// ��ȡʱ��
	double t = GenInitialMatrix[2][0]; //��ǰʱ�������Ϊ��3�е�1��

	// ��ʼ����ǰʱ�̵�abc���ѹ����
	for (int k=0; k<3; k++) {
		Vsabc[k] = GenInitialMatrix[2][ptr+k];
		Vsabc_his[k] = GenInitialMatrix[1][ptr+k];
		Vsabc_his2[k] = GenInitialMatrix[0][ptr+k];
		Vrabc[k] = GenInitialMatrix[2][ptr+3+k];
		Vrabc_his[k] = GenInitialMatrix[1][ptr+3+k];
		Vrabc_his2[k] = GenInitialMatrix[0][ptr+3+k];
		Isabc[k] = GenInitialMatrix[2][ptr+6+k];
		Irabc[k] = GenInitialMatrix[2][ptr+9+k];
	}

	// ��еת�ء�ת�ٺ�ת�ӽ�
	TL = GenInitialMatrix[2][ptr+12];
	wr = GenInitialMatrix[2][ptr+14];
	wr_his = GenInitialMatrix[1][ptr+14];
	 curAngle_s = w*Wbase*t;
	//curAngle_s = 0; // curAngle_s�ĳ�ֵ��������ѡȡ
	double Theta = GenInitialMatrix[2][ptr+13];
	curAngle_r = curAngle_s - Theta;

	// ��ʼ����ѹ������dq����
	double dqMatrix[3][8];

	//��abc�����µĶ��ӵ�ѹ����ת��Ϊ����ֵ
	for (int i=0;i<3;i++)
	{
		for (int j=ptr;j<ptr+3;j++)
			GenInitialMatrix[i][j] /= Vbase;

		for (int j=ptr+3;j<ptr+6;j++)
			GenInitialMatrix[i][j] /= (Vbase/turnratio);

		for (int j=ptr+6;j<ptr+9;j++)
			GenInitialMatrix[i][j] /= Ibase;

		for (int j=ptr+9;j<ptr+12;j++)
			GenInitialMatrix[i][j] /= (Ibase*turnratio);
	}

	// ����dq����
	for (int i=0; i<3; i++)
	{
		// ����park�任�Ƕ�
		double theta_s = w*Wbase*GenInitialMatrix[i][0];
		double theta_r = theta_s - GenInitialMatrix[i][ptr+13];

		// ����park�任����
		parkTransMatrix(Ps,theta_s);
		parkTransMatrix(Pr,theta_r);

		// park�任
		// ���ӵ�ѹ
		for (int j=0; j<2; j++)
		{
			dqMatrix[i][j] = Ps(j+1,1)*GenInitialMatrix[i][ptr] + Ps(j+1,2)*GenInitialMatrix[i][ptr+1]
				+ Ps(j+1,3)*GenInitialMatrix[i][ptr+2];
		}
		// ת�ӵ�ѹ
		for (int j=2; j<4; j++)
		{
			dqMatrix[i][j] = Pr(j-1,1)*GenInitialMatrix[i][ptr+3] + Pr(j-1,2)*GenInitialMatrix[i][ptr+4]
				+ Pr(j-1,3)*GenInitialMatrix[i][ptr+5];
		}
		// ���ӵ���
		for (int j=4; j<6; j++)
		{
			dqMatrix[i][j] = Ps(j-3,1)*GenInitialMatrix[i][ptr+6] + Ps(j-3,2)*GenInitialMatrix[i][ptr+7]
				+ Ps(j-3,3)*GenInitialMatrix[i][ptr+8];
		}
		// ת�ӵ�ѹ
		for (int j=6; j<8; j++)
		{
			dqMatrix[i][j] = Pr(j-5,1)*GenInitialMatrix[i][ptr+9] + Pr(j-5,2)*GenInitialMatrix[i][ptr+10]
				+ Pr(j-5,3)*GenInitialMatrix[i][ptr+11];
		}
	}
	
	// dq������ֵ
	for (int k=0; k<2; k++) {
		Vsdq[k] = dqMatrix[2][k];
		Vsdq_his[k] = dqMatrix[1][k];
		Vsdq_his2[k] = dqMatrix[0][k];
		Isdq[k] = dqMatrix[2][k+4];
		Isdq_his[k] = dqMatrix[1][k+4];
		Isdq_his2[k] = dqMatrix[0][k+4];
		Vrdq[k] = dqMatrix[2][k+2];
		Vrdq_his[k] = dqMatrix[1][k+2];
		Vrdq_his2[k] = dqMatrix[0][k+2];
		Irdq[k] = dqMatrix[2][k+6];
		Irdq_his[k] = dqMatrix[1][k+6];
		Irdq_his2[k] = dqMatrix[0][k+6];
	}

	// �ڵ㵼�ɾ������
	calculateCoefficientMatrix();
	Te_his=0.5*Xm*(Irdq_his2[0]*Isdq_his2[1]-Irdq_his2[1]*Isdq_his2[0]);
	Te=0.5*Xm*(Irdq_his[0]*Isdq_his[1]-Irdq_his[1]*Isdq_his[0]);

	ptr += 15;
}

// �ڷ���������ģʽ�£�������������ת��
void WoundInducGenerator::calculateTwind()
{
#ifdef WIND_VELOCITY_DATA_INPUT
	setWindVelocity();
#endif
	TSR=2.237*vWind/(wr*Wbase/GR);
	Cp=0.5*(TSR-0.022*beta*beta-5.6)*exp(-0.17*TSR);
	Pwind=0.5*Cp*rArea*(vWind*vWind*vWind)*p_air*GE/Sbase;
	Twind=Pwind/wr;
}

//������ת��, By Gao Haixiang
void WoundInducGenerator::saveMachineWr(double** machineWrMatrix, int& ptr, int counter)
{
	machineWrMatrix[counter-1][ptr] = wr;
	ptr++;
}

//	��ÿ��ʱ���ķ������ݱ�����ÿ̨����� ,By Gao Haixiang
void WoundInducGenerator::getWindVelocityData(double** VwMatrix, int rows, int& ptr)
{
	WindVelocityVector = new double[rows];
	for (int i=0;i<rows;i++)
	{
		WindVelocityVector[i] = VwMatrix[i][ptr];
	}
	ptr++;
}

// �趨����, By Gao Haixiang
void WoundInducGenerator::setWindVelocity()
{
	vWind = WindVelocityVector[WindVelocityCounter];
	WindVelocityCounter++;
}