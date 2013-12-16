#include "Motor15.h"

#include <cmath>
#include <iostream>
using namespace std;

Motor15::Motor15(int firstNode,int lastNode)
{
	//�Զ���Ԫ��
	isUserDef = 1;
	need_NEC=1;
	type=15;
	nPort=30;
	nodeNumber=new int[nPort];//chenlj 20090718 added
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}
	//dq0����ϵ�µĵ������
	Rp = 0.1043;
	Rr1 = 0.0188;
	Rrt = 0.0188;
	//dq0����ϵ�µĵ翹����
	Xl1 = 1.2239;
	Xl2 = -0.0163;
	Xlr1 = 0.3101;
	Xm1 = 7.1276;
	Xlt1 = 1.2032;
	Xlt2 = -0.1323;
	Xlrt = 0.2565;
	Xmt = 1.2399;
	Xdqm1 = 0;
	Xdqmt = 0;
	Xs0 = 1.8051;
	//�м����
	Xs1 = Xl1 + Xm1; 
	Xsm1 = Xl2 + Xm1;
	Xr1 = Xlr1 + Xm1;
	Xst = Xlt1 + Xmt;
	Xsmt = Xlt2 + Xmt;
	Xrt = Xlrt + Xmt;

	//���Ӳ��ѹ�������ͽ�Ƶ�ʵĻ�ֵ
	Vbase = 2400;//V,��Чֵ
	Ibase = 684;//A,��Чֵ
	wbase = 2*PI*20;//fbase=20Hz;

	//��е�����
	Je = 25000;//ת������
	p = 6;//������

	//�ɿ˱任�����������
	P.ResizeTo(1,9,1,15);//(9,15)
	invP.ResizeTo(1,15,1,9);//(15,9)

	for (int k=0;k<9;k++)
	{
		for (int i=0;i<5;i++)
		{
			P_1[k][i] = 0;
		}
	}

	for (int k=0;k<15;k++)
	{
		for (int i=0;i<3;i++)
		{
			invP_1[k][i] = 0;
		}
	}

	/*     ���λ��ַ���Ҫʹ�õ�ϵ������      */
	coff = 39.0/41.0;//alpha
	T=1.0/(deltaT*wbase);
	//
	A.ResizeTo(1,11,1,11);//(11,11)
	B.ResizeTo(1,11,1,11);//(11,11)
	AA.ResizeTo(1,11,1,11);//(11,11)
	BB.ResizeTo(1,11,1,11);//(11,11)
	//
	Rss.ResizeTo(1,9,1,9);//AA(1:9,1:9)
	Rsr.ResizeTo(1,9,1,2);//AA(1:9,10:11)
	Rrs.ResizeTo(1,2,1,9);//AA(10:11,1:9)
	Rrr.ResizeTo(1,2,1,2);//AA(10:11,10:11)
	//
	Bss.ResizeTo(1,9,1,9);//BB(1:9,1:9)
	Bsr.ResizeTo(1,9,1,2);//BB(1:9,10:11)
	Brs.ResizeTo(1,2,1,9);//BB(10:11,1:9)
	Brr.ResizeTo(1,2,1,2);//BB(10:11,10:11)

	//����ϵ������ʹ��ת����abcde����ϵʱ�ܽ���
	Rs.ResizeTo(1,9,1,9);//
	Rs1.ResizeTo(1,9,1,9);//
	Rs2.ResizeTo(1,9,1,9);//Rs2=Rs-Rs1;

	//
	Gs1.ResizeTo(1,9,1,9);//
	Gs2.ResizeTo(1,9,1,9);//
	Gs3.ResizeTo(1,9,1,9);//
	Gs4.ResizeTo(1,9,1,2);//

	//
	Gr1.ResizeTo(1,2,1,2);//
	Gr2.ResizeTo(1,2,1,9);//
	Gr3.ResizeTo(1,2,1,9);//
	Gr4.ResizeTo(1,2,1,2);//

	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=0;i<30;i++)
	{
		Vnode[i] = 0;//�ڵ��ѹ
	}
	for (int i=0;i<15;i++)
	{
		Vabcde[i] = 0;//֧·��ѹ
		Iabcde[i] = 0;//֧·����
		Vabcde_1[i] = 0;
		Iabcde_1[i] = 0;
		Vabcde_2[i] = 0;
		Iabcde_2[i] = 0;
	}
	
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=0;i<9;i++)
	{
		Vdqs[i] = 0;
		Idqs[i] = 0;
		Idqs_his[i] = 0;//Ԥ����
		Idqs_his2[i] = 0;//Ԥ����
		Idqs_forcast[i] = 0;//Ԥ��ֵ
	}
	
	//dq0����ϵ�¶��Ӳ�Ĵ���
	psi.ResizeTo(1,11);
	Isr.ResizeTo(1,11);//�������ʱ�õĵ���
	for (int k=0;k<6;k++)
	{
		psi_1[k] = 0;
	}
	//dq0����ϵ��ת�Ӳ�ĵ�ѹ����
	Vr[0] = 0;
	Ir[0] = 0;
	Ir_his[0] = 0;
	Vr[1] = 0;
	Ir[1] = 0;
	Ir_his[1] = 0;
	/*  ��е��������  */
	speed0 = 0;;//ת��ת�ٳ�ʼֵ������ֵ��
	wr0 = speed0*wbase;//ת�ӽ�Ƶ�ʳ�ʼֵ	
	wr = 0;//ת��ת��
	wr_his = 0;//wr����ʷֵ
	theta = 0;//ͬ������ϵd������a��ĽǶȣ�Park�任ʱʹ�ã�
	Te = 0;//���ת��
	Te_his = 0;
	Tm = 0;//���ػ�еת��
	Tm_his = 0;

	nortonEquivalentConductance = 0;
	nortonEquivalentConductance_his = 0;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���
	for (int i=0;i<15;i++)
	{
		nortonEquivalentCurrent[i] = 0;//abcde�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_1[i] = 0;
		nortonEquivalentCurrent_2[i] = 0;
	}
	for (int i=0;i<9;i++)
	{
		nortonEquivalentCurrent_dq0[i] = 0;//dq0�����µ�ŵ�ٵ�Ч����������ֵ
	}

	//�������B
	B(1,1)=Xs1;
	B(1,4)=Xsm1;
	B(1,5)=Xdqm1;
	B(1,7)=Xsm1;
	B(1,8)=Xdqm1;
	B(1,10)=Xm1;

	B(2,2)=Xs1;
	B(2,4)=-Xdqm1;
	B(2,5)=Xsm1;
	B(2,7)=-Xdqm1;
	B(2,8)=Xsm1;
	B(2,11)=Xm1;

	B(3,3)=Xs0;

	B(4,1)=Xsm1;
	B(4,2)=Xdqm1;
	B(4,4)=Xs1;
	B(4,7)=Xsm1;
	B(4,8)=Xdqm1;
	B(4,10)=Xm1;

	B(5,1)=-Xdqm1;
	B(5,2)=Xsm1;
	B(5,5)=Xs1;
	B(5,7)=-Xdqm1;
	B(5,8)=Xsm1;
	B(5,11)=Xm1;

	B(6,6)=Xs0;

	B(7,1)=Xsm1;
	B(7,2)=Xdqm1;
	B(7,4)=Xsm1;
	B(7,5)=Xdqm1;
	B(7,7)=Xs1;
	B(7,10)=Xm1;

	B(8,1)=-Xdqm1;
	B(8,2)=Xsm1;
	B(8,4)=-Xdqm1;
	B(8,5)=Xsm1;
	B(8,8)=Xs1;
	B(8,11)=Xm1;

	B(9,9)=Xs0;
  
	B(10,1)=Xm1;
	B(10,4)=Xm1;
	B(10,7)=Xm1;
	B(10,10)=Xr1;
  
	B(11,2)=Xm1;
	B(11,5)=Xm1;
	B(11,8)=Xm1;
	B(11,11)=Xr1;
  
	//B_1
	B_1[0][0]=B(1,1);
	B_1[0][1]=B(1,4);
	B_1[0][2]=B(1,5);
	B_1[0][3]=B(1,7);
	B_1[0][4]=B(1,8);
	B_1[0][5]=B(1,10);
 
	B_1[1][0]=B(2,2);
	B_1[1][1]=B(2,4);
	B_1[1][2]=B(2,5);
	B_1[1][3]=B(2,7);
	B_1[1][4]=B(2,8);
	B_1[1][5]=B(2,11);
 
	B_1[2][0]=B(4,1);
	B_1[2][1]=B(4,2);
	B_1[2][2]=B(4,4);
	B_1[2][3]=B(4,7);
	B_1[2][4]=B(4,8);
	B_1[2][5]=B(4,10);

	B_1[3][0]=B(5,1);
	B_1[3][1]=B(5,2);
	B_1[3][2]=B(5,5);
	B_1[3][3]=B(5,7);
	B_1[3][4]=B(5,8);
	B_1[3][5]=B(5,11);

	B_1[4][0]=B(7,1);
	B_1[4][1]=B(7,2);
	B_1[4][2]=B(7,4);
	B_1[4][3]=B(7,5);
	B_1[4][4]=B(7,7);
	B_1[4][5]=B(7,10);

	B_1[5][0]=B(8,1);
	B_1[5][1]=B(8,2);
	B_1[5][2]=B(8,4);
	B_1[5][3]=B(8,5);
	B_1[5][4]=B(8,8);
	B_1[5][5]=B(8,11);

	for (int i=0;i<6;i++)
	{
		for (int j=0;j<6;j++)
		{
			B_1[i][j] /= wbase;
		}
	}

	//����Gr1,Gr2,Gr3,Gr4 & nortonEquivalentConductance,Gs1
	Gr1(1,1) = -1/(Rr1+(1+coff)*T*Xr1);
	Gr1(2,2) = Gr1(1,1);

	Gr2(1,1) = -1/(Rr1+(1+coff)*T*Xr1)*(1+coff)*T*Xm1;//Gr1(1,1)*(1+coff)*T*Xm1;
	Gr2(2,2) = Gr2(1,1);
	Gr2(1,4) = Gr2(1,1);
	Gr2(2,5) = Gr2(1,1);
	Gr2(1,7) = Gr2(1,1);
	Gr2(2,8) = Gr2(1,1);

	Gr3(1,1) = 1/(Rr1+(1+coff)*T*Xr1)*(1+coff)*T*Xm1;//-Gr2(1,1)
	Gr3(2,2) = Gr3(1,1);
	Gr3(1,4) = Gr3(1,1);
	Gr3(2,5) = Gr3(1,1);
	Gr3(1,7) = Gr3(1,1);
	Gr3(2,8) = Gr3(1,1);

	Gr4(1,1) = -1/(Rr1+(1+coff)*T*Xr1)*(coff*Rr1-(1+coff)*T*Xr1);
	Gr4(2,2) = Gr4(1,1);

	nortonEquivalentConductance = 1/(Rp+(1+coff)*T*Xs1-((1+coff)*T*Xm1)*((1+coff)*T*Xm1)/(Rr1+(1+coff)*T*Xr1));
	for (int i=1;i<=9;i++)
	{
		Gs1(i,i) = nortonEquivalentConductance;
	}

	temp11 = Xm1*Xm1/(Rr1+(1+coff)*T*Xr1)*(1+coff)*T;
	temp12 = (1+coff)*T*Xsm1;
	temp13 = ((1+coff)*T*Xm1)*((1+coff)*T*Xm1)/(Rr1+(1+coff)*T*Xr1);
	temp14 = (1+coff)*T*Xdqm1;
	temp15 = (1+coff)*T*Xs0;
	temp16 = (1+coff)*T*Xs1;
	temp17 = (1+coff)*T*Xm1;
	temp18 = (1+coff)*T*Xr1;

	//����Gs2,Gs3,Gs4
	calculateCoefficientMatrix();

	for (int i=0;i<9;i++)
	{
		for (int j=0;j<9;j++)
		{
			Gs_1[i][j] = Gs1(i+1,j+1);
			Gs_2[i][j] = Gs2(i+1,j+1);
			Gs_3[i][j] = Gs3(i+1,j+1);
		}
		for (int k=0;k<2;k++)
		{
			Gs_4[i][k] = Gs4(i+1,k+1);
		}
	}

	for (int i=0;i<2;i++)
	{
		for (int j=0;j<9;j++)
		{
			Gr_2[i][j] = Gr2(i+1,j+1);
			Gr_3[i][j] = Gr3(i+1,j+1);
		}
		for (int k=0;k<2;k++)
		{
			Gr_1[i][k] = Gr1(i+1,k+1);
			Gr_4[i][k] = Gr4(i+1,k+1);
		}
	}
}

//��������
Motor15::~Motor15()
{
}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/
//��ʼ��֧·��ѹ����
void Motor15::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	for (int k=0;k<15;k++)
	{
		Vabcde[k] = 0;
		Iabcde[k] = 0;
	}

	Iabcde[0]=1000.0;
	Iabcde[1]=1000.0;
	Iabcde[2]=0.0;
	Iabcde[3]=-1000.0;
	Iabcde[4]=-1000.0;

	//����Park�任����
	double Angle;
	Angle = theta;
	parkTransMatrix(P,Angle);

	//����dq0����ϵ�µĶ��ӵ�ѹ
	for (int i=0;i<9;i++)
	{
		Vdqs[i] = 0;
		for (int j=0;j<15;j++)
		{
			Vdqs[i] += P(i+1,j+1)*Vabcde[j];
		}
	}

	//����dq0����ϵ�µĶ��ӵ���
	for (int i=0;i<9;i++)
	{
		Idqs[i] = 0;
		for (int j=0;j<15;j++)
		{
			Idqs[i] += P(i+1,j+1)*Iabcde[j];
		}
	}

	for (int k=0;k<9;k++)
	{
		Idqs_his2[k] = Idqs[k];
		Idqs_his[k] = Idqs[k];
	}

	//Ԥ����һʱ�̶��ӵ���
	for (int i=0;i<9;i++)
	{
		Idqs_forcast[i] = (202.0/120.0)*Idqs[i]-(44.0/120.0)*Idqs_his[i]-(38.0/120.0)*Idqs_his2[i];
	}

	//����ת�ӵ���
	Ir[0] = 0;
	Ir[1] = 0;

	Ir_his[0] = Ir[0];
	Ir_his[1] = Ir[1];

	for (int i=0;i<2;i++)
	{
		Ir[i] = 0;
		for (int j=0;j<9;j++)
		{
			Ir[i] = Ir[i] + Gr2(i+1,j+1)*Idqs[j] + Gr3(i+1,j+1)*Idqs_his[j];
		}
		Ir[i] = Ir[i] + Gr4(i+1,1)*Ir_his[0] + Gr4(i+1,2)*Ir_his[1];
	}

	//����Inew
	TVectorD tem1,Idqs_forcast_V,Idqs_V,Vdqs_V,Ir_V,Inew;
	tem1.ResizeTo(1,9);
	Idqs_forcast_V.ResizeTo(1,9);
	Idqs_V.ResizeTo(1,9);
	Vdqs_V.ResizeTo(1,9);
	Ir_V.ResizeTo(1,2);
	Inew.ResizeTo(1,9);

	for (int i=1;i<=9;i++)
	{
		Idqs_forcast_V(i) = Idqs_forcast[i-1];
		Idqs_V(i) = Idqs[i-1];
		Vdqs_V(i) = Vdqs[i-1];
	}
	Ir_V(1) = Ir[0];
	Ir_V(2) = Ir[1];

	tem1 = Gs2*Idqs_forcast_V + Gs3*Idqs_V + Gs4*Ir_V + coff*Vdqs_V;
	Inew = Gs1*tem1;

	//����abc�����µ�ŵ�ٵ�ֵ����
	for (int i=0;i<9;i++)
	{
		nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	}

	Angle = wr*deltaT + theta;
	invParkTransMatrix(invP,Angle);

	for (int i=0;i<15;i++)
	{
		nortonEquivalentCurrent[i] = 0;
		for (int j=0;j<9;j++)
		{
			nortonEquivalentCurrent[i] += invP(i+1,j+1)*nortonEquivalentCurrent_dq0[j];
		}
	}
	
	//���¼���Idqs��Idqs_his��Idqs_his2
	for (int k=0;k<9;k++)
	{
		Idqs[k] = Inew(k+1);
		Idqs_his[k] = Idqs[k];
		Idqs_his2[k] = Idqs[k];
	}

	//����е����
	//����Isr
	for (int i=1;i<=9;i++)
	{
		Isr(i) = Idqs[i-1];
	}
	for (int i=1;i<=2;i++)
	{
		Isr(i+9) = Ir[i-1];
	}

	//�������psi
	psi = B*Isr;

	for (int i=1;i<=11;i++)
	{
		psi(i) /= wbase;
	}

	//������ת��Te
	Te = (psi(1)*Idqs[1]-psi(2)*Idqs[0])*15.0/2*p/1e6;

	//����ת�ӽ�Ƶ��wr
	wr_his = wr;
	wr = wr_his + (Te-Tm)*deltaT/Je;

	//����ת�ӽ�theta
	theta = theta+(wr-wr0)*deltaT;

	//��PSCAD���ݳ�ʼ�����ӵ�ѹ����
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<15;i++)
	{
		Iabcde[i] = initialCurrentArray[ptr+i];
	}
	ptr+=15;
}

//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
void Motor15::readNodeVoltage(TVectorD& nodeVoltageArray)
{//�ӽڵ��ѹ�����ж��綯���˽ڵ��ѹ
	for (int i=0;i<30;i++)
	{
		if (nodeNumber[i] == 0)
		{
			Vnode[i] = 0;
		}else
		{
			Vnode[i] = nodeVoltageArray(nodeNumber[i]);
		}
	}
}

//����֧·��ѹ
void Motor15::calculateBranchVoltage()
{
	for (int i=0;i<15;i++)
	{
		Vabcde_1[i] = Vabcde[i];
		Vabcde[i] = Vnode[2*i] - Vnode[2*i+1];
	}
}

//����֧·����
void Motor15::calculateBranchCurrent()
{
	for (int i=0;i<15;i++)
	{
		Iabcde_1[i] = Iabcde[i];
		Iabcde[i] = nortonEquivalentConductance*Vabcde[i] + nortonEquivalentCurrent[i];
	}
}

//����֧·��ŵ�ٵ�Ч��·�еĵ�����
void Motor15::calculateNortonEquivalentCurrent(double time)
{
	//calculateCoefficientMatrix();
	calculateCoefficientMatrix_1();
	calculateDq0Results(time);

	//TVectorD tem1,Idqs_forcast_V,Idqs_V,Vdqs_V,Ir_V,Inew;
	//tem1.ResizeTo(1,9);
	//Idqs_forcast_V.ResizeTo(1,9);
	//Idqs_V.ResizeTo(1,9);
	//Vdqs_V.ResizeTo(1,9);
	//Ir_V.ResizeTo(1,2);
	//Inew.ResizeTo(1,9);

	//for (int i=1;i<=9;i++)
	//{
	//	Idqs_forcast_V(i) = Idqs_forcast[i-1];
	//	Idqs_V(i) = Idqs[i-1];
	//	Vdqs_V(i) = Vdqs[i-1];
	//}
	//Ir_V(1) = Ir[0];
	//Ir_V(2) = Ir[1];

	//tem1 = Gs2*Idqs_forcast_V + Gs3*Idqs_V + Gs4*Ir_V + coff*Vdqs_V;
	//Inew = Gs1*tem1;

	//for (int i=0;i<9;i++)
	//{
	//	nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	//}

	double tem1[9];
	for (int i=0;i<9;i++)
	{
		tem1[i] = 0;
		for (int j=0;j<9;j++)
		{
			tem1[i] = tem1[i] + Gs_2[i][j]*Idqs_forcast[j] + Gs_3[i][j]*Idqs[j];
		}
		tem1[i] = tem1[i] + Gs_4[i][0]*Ir[0] + Gs_4[i][1]*Ir[1] + coff*Vdqs[i];
	}
	for (int i=0;i<9;i++)
	{
		nortonEquivalentCurrent_dq0[i] = 0;
		for (int j=0;j<9;j++)
		{
			nortonEquivalentCurrent_dq0[i] += Gs_1[i][j]*tem1[j];
		}
	}

	double Angle;
	Angle = wr0*time + wr*deltaT + theta;
	//invParkTransMatrix(invP,Angle);
	invParkTransMatrix(invP_1,Angle);

	for (int k=0;k<15;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//for (int i=0;i<15;i++)
	//{
	//	nortonEquivalentCurrent[i] = 0;
	//	for (int j=0;j<9;j++)
	//	{
	//		nortonEquivalentCurrent[i] += invP(i+1,j+1)*nortonEquivalentCurrent_dq0[j];
	//	}
	//}

	for (int i=0;i<3;i++)
	{
		for (int j=0;j<5;j++)
		{
			nortonEquivalentCurrent[5*i+j] = invP_1[5*i+j][0]*nortonEquivalentCurrent_dq0[3*i]+invP_1[5*i+j][1]*nortonEquivalentCurrent_dq0[3*i+1]+invP_1[5*i+j][2]*nortonEquivalentCurrent_dq0[3*i+2];
		}
	}

	calculateOmega();
}

//�γɽڵ�ŵ�ٵ�Ч��������
void Motor15::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{
	int k,m;
	for (int i=0;i<15;i++)
	{
		k = nodeNumber[2*i];
		m = nodeNumber[2*i+1];

		nodeNortonEquivalentCurrentArray(k) -= nortonEquivalentCurrent[i];
		nodeNortonEquivalentCurrentArray(m) += nortonEquivalentCurrent[i];


		//if (isSeries(k,m)==1)
		//{
		//	nodeNortonEquivalentCurrentArray(k) -= nortonEquivalentCurrent[i];
		//	nodeNortonEquivalentCurrentArray(m) += nortonEquivalentCurrent[i];
		//}
		//else
		//{
		//	if (m==0)
		//	{
		//		nodeNortonEquivalentCurrentArray(k) -= nortonEquivalentCurrent[i];
		//	} 
		//	else
		//	{
		//		nodeNortonEquivalentCurrentArray(m) += nortonEquivalentCurrent[i];
		//	}
		//}
	}
}
void Motor15::formConductanceMatrix(TMatrixD &conductanceMatrix)
{
	int k,m;
	for (int i=0;i<15;i++)
	{
		k = nodeNumber[2*i];
		m = nodeNumber[2*i+1];

		conductanceMatrix(k,k) += nortonEquivalentConductance;
		conductanceMatrix(m,m) += nortonEquivalentConductance;
		conductanceMatrix(k,m) -= nortonEquivalentConductance;
		conductanceMatrix(m,k) -= nortonEquivalentConductance;


		//if (isSeries(k,m)==1)
		//{
		//	conductanceMatrix(k,k) += nortonEquivalentConductance;
		//	conductanceMatrix(m,m) += nortonEquivalentConductance;
		//	conductanceMatrix(k,m) -= nortonEquivalentConductance;
		//	conductanceMatrix(m,k) -= nortonEquivalentConductance;
		//} 
		//else
		//{
		//	if (m==0)
		//	{
		//		conductanceMatrix(k,k) += nortonEquivalentConductance;
		//	} 
		//	else
		//	{
		//		conductanceMatrix(m,m) += nortonEquivalentConductance;
		//	}
		//}
	}
}

//����֧·����
void Motor15::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{
	for (int i=0;i<15;i++)
	{
		branchCurrentMatrix(counter,ptr+i) = Iabcde[i];
	}
	ptr+=15;
}

//����֧·����
void Motor15::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int i=0;i<15;i++)
	{
		branchCurrentMatrix_1[counter-1][ptr+i] = Iabcde[i];
	}
	ptr+=15;
}

//������ֵ����֧·�ĵ�ѹ�������в�ֵ
void Motor15::interpolate(double ratio)
{
	//��ѹ����������ֵ
	for (int k=0;k<15;k++)
	{
		Vabcde[k] = (1-ratio)*Vabcde_1[k] + ratio*Vabcde[k];
		Iabcde[k] = (1-ratio)*Iabcde_1[k] + ratio*Iabcde[k];
	}

	//ŵ�ٵ�ֵ������ֵ
	double temp1,temp2;
	for (int k=0;k<15;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void Motor15::updateResult(int updateTypeNum)
{
	double V_temp[15],I_temp[15];//��ʱ��������������ʱʹ��

	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for (int i=0;i<15;i++)
		{
			Vabcde_2[i] = Vabcde_1[i];
			Iabcde_2[i] = Iabcde_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://��_2��������ֵ����_1������
		for (int i=0;i<15;i++)
		{
			Vabcde_1[i] = Vabcde_2[i];
			Iabcde_1[i] = Iabcde_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://ŵ�ٵ�ֵ�����洢����_1��������ֵ����_2������
		for (int i=0;i<15;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"��������ȷ�ĸ������ͱ�ţ�"<<endl;
		exit(1);
	}
}

/************************************************************************/
/*                �綯���ڲ����õĺ���                                  */
/************************************************************************/
//����ϵ������AA,BB,Rss...Brr,Reqs...inv_Rr,Gs1...Gs4,Gr1...Gr4,...
void Motor15::calculateCoefficientMatrix()
{
	double w;//ת�ӽ�Ƶ�ʱ���ֵ
	w = wr/wbase;

	//����Gs2
	Gs2(1,2) = w*Xs1-w*temp11;
	Gs2(1,4) = -w*Xdqm1-temp12+temp13;
	Gs2(1,5) = w*Xsm1-temp14-w*temp11;
	Gs2(1,7) = Gs2(1,4);
	Gs2(1,8) = Gs2(1,5);

	Gs2(2,1) = -Gs2(1,2);
	Gs2(2,4) = -Gs2(1,5);
	Gs2(2,5) = Gs2(1,4);
	Gs2(2,7) = Gs2(2,4);
	Gs2(2,8) = Gs2(2,5);

	Gs2(3,3) = -temp15+temp16-temp13;

	Gs2(4,1) = Gs2(1,4);
	Gs2(4,2) = Gs2(1,5);
	Gs2(4,5) = Gs2(1,2);
	Gs2(4,7) = Gs2(1,4);
	Gs2(4,8) = Gs2(1,5);

	Gs2(5,1) = Gs2(2,4);
	Gs2(5,2) = Gs2(2,5);
	Gs2(5,4) = Gs2(2,1);
	Gs2(5,7) = Gs2(2,4);
	Gs2(5,8) = Gs2(2,5);

	Gs2(6,6) = Gs2(3,3);

	Gs2(7,1) = Gs2(1,4);
	Gs2(7,2) = Gs2(1,5);
	Gs2(7,4) = Gs2(1,4);
	Gs2(7,5) = Gs2(1,5);
	Gs2(7,8) = Gs2(1,2);

	Gs2(8,1) = Gs2(2,4);
	Gs2(8,2) = Gs2(2,5);
	Gs2(8,4) = Gs2(2,4);
	Gs2(8,5) = Gs2(2,5);
	Gs2(8,7) = Gs2(2,1);

	Gs2(9,9) = Gs2(3,3);

	//����Gs3
	Gs3(1,1) = -coff*Rp+temp16-temp13;
	Gs3(1,2) = coff*w*Xs1+w*temp11;
	Gs3(1,4) = -coff*w*Xdqm1+temp12-temp13;
	Gs3(1,5) = coff*w*Xsm1+temp14+w*temp11;
	Gs3(1,7) = Gs3(1,4);
	Gs3(1,8) = Gs3(1,5);

	Gs3(2,1) = -Gs3(1,2);
	Gs3(2,2) = Gs3(1,1);
	Gs3(2,4) = -Gs3(1,5);
	Gs3(2,5) = Gs3(1,4);
	Gs3(2,7) = Gs3(2,4);
	Gs3(2,8) = Gs3(2,5);

	Gs3(3,3) = -coff*Rp+temp15;

	Gs3(4,1) = Gs3(1,4);
	Gs3(4,2) = Gs3(1,5);
	Gs3(4,4) = Gs3(1,1);
	Gs3(4,5) = Gs3(1,2);
	Gs3(4,7) = Gs3(1,4);
	Gs3(4,8) = Gs3(1,5);

	Gs3(5,1) = Gs3(2,4);
	Gs3(5,2) = Gs3(2,5);
	Gs3(5,4) = Gs3(2,1);
	Gs3(5,5) = Gs3(2,2);
	Gs3(5,7) = Gs3(2,4);
	Gs3(5,8) = Gs3(2,5);

	Gs3(6,6) = Gs3(3,3);

	Gs3(7,1) = Gs3(1,4);
	Gs3(7,2) = Gs3(1,5);
	Gs3(7,4) = Gs3(1,4);
	Gs3(7,5) = Gs3(1,5);
	Gs3(7,7) = Gs3(1,1);
	Gs3(7,8) = Gs3(1,2);

	Gs3(8,1) = Gs3(2,4);
	Gs3(8,2) = Gs3(2,5);
	Gs3(8,4) = Gs3(2,4);
	Gs3(8,5) = Gs3(2,5);
	Gs3(8,7) = Gs3(2,1);
	Gs3(8,8) = Gs3(2,2);

	Gs3(9,9) = Gs3(3,3);

	//����Gs4
	Gs4(1,1) = temp17+temp17/(Rr1+temp18)*(coff*Rr1-temp18);
	Gs4(1,2) = coff*w*Xm1-w*Xm1/(Rr1+temp18)*(coff*Rr1-temp18);

	Gs4(2,1) = (-1.0)*Gs4(1,2);
	Gs4(2,2) = Gs4(1,1);

	Gs4(4,1) = Gs4(1,1);
	Gs4(4,2) = Gs4(1,2);

	Gs4(5,1) = Gs4(2,1);
	Gs4(5,2) = Gs4(2,2);

	Gs4(7,1) = Gs4(1,1);
	Gs4(7,2) = Gs4(1,2);

	Gs4(8,1) = Gs4(2,1);
	Gs4(8,2) = Gs4(2,2);
}

void Motor15::calculateCoefficientMatrix_1()
{
	double w;//ת�ӽ�Ƶ�ʱ���ֵ
	w = wr/wbase;

	//����Gs2
	Gs_2[0][1] = w*Xs1-w*temp11;
	Gs_2[0][3] = -w*Xdqm1-temp12+temp13;
	Gs_2[0][4] = w*Xsm1-temp14-w*temp11;
	Gs_2[0][6] = Gs_2[0][3];
	Gs_2[0][7] = Gs_2[0][4];

	Gs_2[1][0] = -Gs_2[0][1];
	Gs_2[1][3] = -Gs_2[0][4];
	Gs_2[1][4] = Gs_2[0][3];
	Gs_2[1][6] = Gs_2[1][3];
	Gs_2[1][7] = Gs_2[1][4];

	Gs_2[2][2] = -temp15+temp16-temp13;

	Gs_2[3][0] = Gs_2[0][3];
	Gs_2[3][1] = Gs_2[0][4];
	Gs_2[3][4] = Gs_2[0][1];
	Gs_2[3][6] = Gs_2[0][3];
	Gs_2[3][7] = Gs_2[0][4];

	Gs_2[4][0] = Gs_2[1][3];
	Gs_2[4][1] = Gs_2[1][4];
	Gs_2[4][3] = Gs_2[1][0];
	Gs_2[4][6] = Gs_2[1][3];
	Gs_2[4][7] = Gs_2[1][4];

	Gs_2[5][5] = Gs_2[2][2];

	Gs_2[6][0] = Gs_2[0][3];
	Gs_2[6][1] = Gs_2[0][4];
	Gs_2[6][3] = Gs_2[0][3];
	Gs_2[6][4] = Gs_2[0][4];
	Gs_2[6][7] = Gs_2[0][1];

	Gs_2[7][0] = Gs_2[1][3];
	Gs_2[7][1] = Gs_2[1][4];
	Gs_2[7][3] = Gs_2[1][3];
	Gs_2[7][4] = Gs_2[1][4];
	Gs_2[7][6] = Gs_2[1][0];

	Gs_2[8][8] = Gs_2[2][2];

	//����Gs3
	Gs_3[0][0] = -coff*Rp+temp16-temp13;
	Gs_3[0][1] = coff*w*Xs1+w*temp11;
	Gs_3[0][3] = -coff*w*Xdqm1+temp12-temp13;
	Gs_3[0][4] = coff*w*Xsm1+temp14+w*temp11;
	Gs_3[0][6] = Gs_3[0][3];
	Gs_3[0][7] = Gs_3[0][4];

	Gs_3[1][0] = -Gs_3[0][1];
	Gs_3[1][1] = Gs_3[0][0];
	Gs_3[1][3] = -Gs_3[0][4];
	Gs_3[1][4] = Gs_3[0][3];
	Gs_3[1][6] = Gs_3[1][3];
	Gs_3[1][7] = Gs_3[1][4];

	Gs_3[2][2] = -coff*Rp+temp15;

	Gs_3[3][0] = Gs_3[0][3];
	Gs_3[3][1] = Gs_3[0][4];
	Gs_3[3][3] = Gs_3[0][0];
	Gs_3[3][4] = Gs_3[0][1];
	Gs_3[3][6] = Gs_3[0][3];
	Gs_3[3][7] = Gs_3[0][4];

	Gs_3[4][0] = Gs_3[1][3];
	Gs_3[4][1] = Gs_3[1][4];
	Gs_3[4][3] = Gs_3[1][0];
	Gs_3[4][4] = Gs_3[1][1];
	Gs_3[4][6] = Gs_3[1][3];
	Gs_3[4][7] = Gs_3[1][4];

	Gs_3[5][5] = Gs_3[2][2];

	Gs_3[6][0] = Gs_3[0][3];
	Gs_3[6][1] = Gs_3[0][4];
	Gs_3[6][3] = Gs_3[0][3];
	Gs_3[6][4] = Gs_3[0][4];
	Gs_3[6][6] = Gs_3[0][0];
	Gs_3[6][7] = Gs_3[0][1];

	Gs_3[7][0] = Gs_3[1][3];
	Gs_3[7][1] = Gs_3[1][4];
	Gs_3[7][3] = Gs_3[1][3];
	Gs_3[7][4] = Gs_3[1][4];
	Gs_3[7][6] = Gs_3[1][0];
	Gs_3[7][7] = Gs_3[1][1];

	Gs_3[8][8] = Gs_3[2][2];

	//����Gs4
	Gs_4[0][0] = temp17+temp17/(Rr1+temp18)*(coff*Rr1-temp18);
	Gs_4[0][1] = coff*w*Xm1-w*Xm1/(Rr1+temp18)*(coff*Rr1-temp18);

	Gs_4[1][0] = (-1.0)*Gs_4[0][1];
	Gs_4[1][1] = Gs_4[0][0];

	Gs_4[3][0] = Gs_4[0][0];
	Gs_4[3][1] = Gs_4[0][1];

	Gs_4[4][0] = Gs_4[1][0];
	Gs_4[4][1] = Gs_4[1][1];

	Gs_4[6][0] = Gs_4[0][0];
	Gs_4[6][1] = Gs_4[0][1];

	Gs_4[7][0] = Gs_4[1][0];
	Gs_4[7][1] = Gs_4[1][1];
}


//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
void Motor15::calculateDq0Results(double time)
{
	//����Park�任����
	double Angle;
	Angle = wr0*time+theta;
	//parkTransMatrix(P,Angle);
	parkTransMatrix(P_1,Angle);

	//����dq0����ϵ�µĶ��ӵ�ѹ
	//for (int i=0;i<9;i++)
	//{
	//	Vdqs[i] = 0;
	//	for (int j=0;j<15;j++)
	//	{
	//		Vdqs[i] += P(i+1,j+1)*Vabcde[j];
	//	}
	//}
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			Vdqs[3*i+j] = P_1[3*i+j][0]*Vabcde[5*i]+P_1[3*i+j][1]*Vabcde[5*i+1]+P_1[3*i+j][2]*Vabcde[5*i+2]+P_1[3*i+j][3]*Vabcde[5*i+3]+P_1[3*i+j][4]*Vabcde[5*i+4];
		}
	}

	//����dq0����ϵ�µĶ��ӵ���
	//for (int i=0;i<9;i++)
	//{
	//	Idqs_his2[i] = Idqs_his[i];
	//	Idqs_his[i] = Idqs[i];
	//	Idqs[i] = 0;
	//	for (int j=0;j<15;j++)
	//	{
	//		Idqs[i] += P(i+1,j+1)*Iabcde[j];
	//	}
	//}

	for (int i=0;i<9;i++)
	{
		Idqs_his2[i] = Idqs_his[i];
		Idqs_his[i] = Idqs[i];
	}
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			Idqs[3*i+j] = P_1[3*i+j][0]*Iabcde[5*i]+P_1[3*i+j][1]*Iabcde[5*i+1]+P_1[3*i+j][2]*Iabcde[5*i+2]+P_1[3*i+j][3]*Iabcde[5*i+3]+P_1[3*i+j][4]*Iabcde[5*i+4];
		}
	}

	//Ԥ����һʱ�̶��ӵ���
	for (int i=0;i<9;i++)
	{
		Idqs_forcast[i] = (140.0/120.0)*Idqs[i]-(10.0/120.0)*Idqs_his[i]-(10.0/120.0)*Idqs_his2[i];
	}

	//����ת�ӵ���
	Ir_his[0] = Ir[0];
	Ir_his[1] = Ir[1];

	//for (int i=0;i<2;i++)
	//{
	//	Ir[i] = 0;
	//	for (int j=0;j<9;j++)
	//	{
	//		Ir[i] = Ir[i] + Gr2(i+1,j+1)*Idqs[j] + Gr3(i+1,j+1)*Idqs_his[j];
	//	}
	//	Ir[i] = Ir[i] + Gr4(i+1,1)*Ir_his[0] + Gr4(i+1,2)*Ir_his[1];
	//}

	for (int i=0;i<2;i++)
	{
		Ir[i] = 0;
		for (int j=0;j<9;j++)
		{
			Ir[i] = Ir[i] + Gr_2[i][j]*Idqs[j] + Gr_3[i][j]*Idqs_his[j];
		}
		Ir[i] = Ir[i] + Gr_4[i][0]*Ir_his[0] + Gr_4[i][1]*Ir_his[1];
	}
}

//����е���̣�����ת��ת��
void Motor15::calculateOmega()
{
	////����Isr
	//for (int i=1;i<=9;i++)
	//{
	//	Isr(i) = Idqs[i-1];
	//}
	//for (int i=1;i<=2;i++)
	//{
	//	Isr(i+9) = Ir[i-1];
	//}

	////�������psi
	//psi = B*Isr;

	//for (int i=1;i<=11;i++)
	//{
	//	psi(i) /= wbase;
	//}

	psi_1[0] = B_1[0][0]*Idqs[0] + B_1[0][1]*Idqs[3] + B_1[0][2]*Idqs[4] + B_1[0][3]*Idqs[6] + B_1[0][4]*Idqs[7] + B_1[0][5]*Ir[0];
	psi_1[1] = B_1[1][0]*Idqs[1] + B_1[1][1]*Idqs[3] + B_1[1][2]*Idqs[4] + B_1[1][3]*Idqs[6] + B_1[1][4]*Idqs[7] + B_1[1][5]*Ir[1];
	psi_1[2] = B_1[2][0]*Idqs[0] + B_1[2][1]*Idqs[1] + B_1[2][2]*Idqs[3] + B_1[2][3]*Idqs[6] + B_1[2][4]*Idqs[7] + B_1[2][5]*Ir[0];
	psi_1[3] = B_1[3][0]*Idqs[0] + B_1[3][1]*Idqs[1] + B_1[3][2]*Idqs[4] + B_1[3][3]*Idqs[6] + B_1[3][4]*Idqs[7] + B_1[3][5]*Ir[1];
	psi_1[4] = B_1[4][0]*Idqs[0] + B_1[4][1]*Idqs[1] + B_1[4][2]*Idqs[3] + B_1[4][3]*Idqs[4] + B_1[4][4]*Idqs[6] + B_1[4][5]*Ir[0];
	psi_1[5] = B_1[5][0]*Idqs[0] + B_1[5][1]*Idqs[1] + B_1[5][2]*Idqs[3] + B_1[5][3]*Idqs[4] + B_1[5][4]*Idqs[7] + B_1[5][5]*Ir[1];

	//������ת��Te
	//Te = (psi(1)*Idqs[1]-psi(2)*Idqs[0]+psi(4)*Idqs[4]-psi(5)*Idqs[3]+psi(7)*Idqs[7]-psi(8)*Idqs[6])*5.0/2*p/1e6;
	Te = (psi_1[0]*Idqs[1]-psi_1[1]*Idqs[0]+psi_1[2]*Idqs[4]-psi_1[3]*Idqs[3]+psi_1[4]*Idqs[7]-psi_1[5]*Idqs[6])*5.0/2*p/1e6;

	//����ת�ӽ�Ƶ��wr
	wr_his = wr;
	wr = wr_his + p*(Te-Tm)*deltaT/Je;

	//����ת�ӽ�theta
	theta = theta+(wr-wr0)*deltaT;
}

//�����ɿ˱任����
void Motor15::parkTransMatrix(TMatrixD &P,double Angle)
{
	P(1,1)=2.0/5.0*cos(Angle);
	P(1,2)=2.0/5.0*cos(Angle-2.0/5.0*PI);
	P(1,3)=2.0/5.0*cos(Angle-4.0/5.0*PI);
	P(1,4)=2.0/5.0*cos(Angle+4.0/5.0*PI);
	P(1,5)=2.0/5.0*cos(Angle+2.0/5.0*PI);

	P(2,1)=-2.0/5.0*sin(Angle);
	P(2,2)=-2.0/5.0*sin(Angle-2.0/5.0*PI);
	P(2,3)=-2.0/5.0*sin(Angle-4.0/5.0*PI);
	P(2,4)=-2.0/5.0*sin(Angle+4.0/5.0*PI);
	P(2,5)=-2.0/5.0*sin(Angle+2.0/5.0*PI);

	P(3,1)=1.0/5.0;
	P(3,2)=1.0/5.0;
	P(3,3)=1.0/5.0;
	P(3,4)=1.0/5.0;
	P(3,5)=1.0/5.0;

	P(4,6)=2.0/5.0*cos(Angle-1.0/15.0*PI);
	P(4,7)=2.0/5.0*cos(Angle-2.0/5.0*PI-1.0/15.0*PI);
	P(4,8)=2.0/5.0*cos(Angle-4.0/5.0*PI-1.0/15.0*PI);
	P(4,9)=2.0/5.0*cos(Angle+4.0/5.0*PI-1.0/15.0*PI);
	P(4,10)=2.0/5.0*cos(Angle+2.0/5.0*PI-1.0/15.0*PI);

	P(5,6)=-2.0/5.0*sin(Angle-1.0/15.0*PI);
	P(5,7)=-2.0/5.0*sin(Angle-2.0/5.0*PI-1.0/15.0*PI);
	P(5,8)=-2.0/5.0*sin(Angle-4.0/5.0*PI-1.0/15.0*PI);
	P(5,9)=-2.0/5.0*sin(Angle+4.0/5.0*PI-1.0/15.0*PI);
	P(5,10)=-2.0/5.0*sin(Angle+2.0/5.0*PI-1.0/15.0*PI);

	P(6,6)=1.0/5.0;
	P(6,7)=1.0/5.0;
	P(6,8)=1.0/5.0;
	P(6,9)=1.0/5.0;
	P(6,10)=1.0/5.0;

	P(7,11)=2.0/5.0*cos(Angle-2.0/15.0*PI);
	P(7,12)=2.0/5.0*cos(Angle-2.0/5.0*PI-2.0/15.0*PI);
	P(7,13)=2.0/5.0*cos(Angle-4.0/5.0*PI-2.0/15.0*PI);
	P(7,14)=2.0/5.0*cos(Angle+4.0/5.0*PI-2.0/15.0*PI);
	P(7,15)=2.0/5.0*cos(Angle+2.0/5.0*PI-2.0/15.0*PI);

	P(8,11)=-2.0/5.0*sin(Angle-2.0/15.0*PI);
	P(8,12)=-2.0/5.0*sin(Angle-2.0/5.0*PI-2.0/15.0*PI);
	P(8,13)=-2.0/5.0*sin(Angle-4.0/5.0*PI-2.0/15.0*PI);
	P(8,14)=-2.0/5.0*sin(Angle+4.0/5.0*PI-2.0/15.0*PI);
	P(8,15)=-2.0/5.0*sin(Angle+2.0/5.0*PI-2.0/15.0*PI);

	P(9,11)=1.0/5.0;
	P(9,12)=1.0/5.0;
	P(9,13)=1.0/5.0;
	P(9,14)=1.0/5.0;
	P(9,15)=1.0/5.0;
}

//�����ɿ˷��任����
void Motor15::invParkTransMatrix(TMatrixD &INRP,double Angle)
{
	INRP(1,1)=cos(Angle);
	INRP(2,1)=cos(Angle-2.0/5.0*PI);
	INRP(3,1)=cos(Angle-4.0/5.0*PI);
	INRP(4,1)=cos(Angle+4.0/5.0*PI);
	INRP(5,1)=cos(Angle+2.0/5.0*PI);

	INRP(1,2)=-sin(Angle);
	INRP(2,2)=-sin(Angle-2.0/5.0*PI);
	INRP(3,2)=-sin(Angle-4.0/5.0*PI);
	INRP(4,2)=-sin(Angle+4.0/5.0*PI);
	INRP(5,2)=-sin(Angle+2.0/5.0*PI);

	INRP(1,3)=1.0;
	INRP(2,3)=1.0;
	INRP(3,3)=1.0;
	INRP(4,3)=1.0;
	INRP(5,3)=1.0;

	INRP(6,4)=cos(Angle-1.0/15.0*PI);
	INRP(7,4)=cos(Angle-2.0/5.0*PI-1.0/15.0*PI);
	INRP(8,4)=cos(Angle-4.0/5.0*PI-1.0/15.0*PI);
	INRP(9,4)=cos(Angle+4.0/5.0*PI-1.0/15.0*PI);
	INRP(10,4)=cos(Angle+2.0/5.0*PI-1.0/15.0*PI);

	INRP(6,5)=-sin(Angle-1.0/15.0*PI);
	INRP(7,5)=-sin(Angle-2.0/5.0*PI-1.0/15.0*PI);
	INRP(8,5)=-sin(Angle-4.0/5.0*PI-1.0/15.0*PI);
	INRP(9,5)=-sin(Angle+4.0/5.0*PI-1.0/15.0*PI);
	INRP(10,5)=-sin(Angle+2.0/5.0*PI-1.0/15.0*PI);

	INRP(6,6)=1.0;
	INRP(7,6)=1.0;
	INRP(8,6)=1.0;
	INRP(9,6)=1.0;
	INRP(10,6)=1.0;

	INRP(11,7)=cos(Angle-2.0/15.0*PI);
	INRP(12,7)=cos(Angle-2.0/5.0*PI-2.0/15.0*PI);
	INRP(13,7)=cos(Angle-4.0/5.0*PI-2.0/15.0*PI);
	INRP(14,7)=cos(Angle+4.0/5.0*PI-2.0/15.0*PI);
	INRP(15,7)=cos(Angle+2.0/5.0*PI-2.0/15.0*PI);

	INRP(11,8)=-sin(Angle-2.0/15.0*PI);
	INRP(12,8)=-sin(Angle-2.0/5.0*PI-2.0/15.0*PI);
	INRP(13,8)=-sin(Angle-4.0/5.0*PI-2.0/15.0*PI);
	INRP(14,8)=-sin(Angle+4.0/5.0*PI-2.0/15.0*PI);
	INRP(15,8)=-sin(Angle+2.0/5.0*PI-2.0/15.0*PI);

	INRP(11,9)=1.0;
	INRP(12,9)=1.0;
	INRP(13,9)=1.0;
	INRP(14,9)=1.0;
	INRP(15,9)=1.0;
}

void Motor15::parkTransMatrix(double P[9][5],double Angle)
{
	P[0][0]=2.0/5.0*cos(Angle);
	P[0][1]=2.0/5.0*cos(Angle-2.0/5.0*PI);
	P[0][2]=2.0/5.0*cos(Angle-4.0/5.0*PI);
	P[0][3]=2.0/5.0*cos(Angle+4.0/5.0*PI);
	P[0][4]=2.0/5.0*cos(Angle+2.0/5.0*PI);

	P[1][0]=-2.0/5.0*sin(Angle);
	P[1][1]=-2.0/5.0*sin(Angle-2.0/5.0*PI);
	P[1][2]=-2.0/5.0*sin(Angle-4.0/5.0*PI);
	P[1][3]=-2.0/5.0*sin(Angle+4.0/5.0*PI);
	P[1][4]=-2.0/5.0*sin(Angle+2.0/5.0*PI);

	P[2][0]=1.0/5.0;
	P[2][1]=1.0/5.0;
	P[2][2]=1.0/5.0;
	P[2][3]=1.0/5.0;
	P[2][4]=1.0/5.0;

	P[3][0]=2.0/5.0*cos(Angle-1.0/15.0*PI);
	P[3][1]=2.0/5.0*cos(Angle-2.0/5.0*PI-1.0/15.0*PI);
	P[3][2]=2.0/5.0*cos(Angle-4.0/5.0*PI-1.0/15.0*PI);
	P[3][3]=2.0/5.0*cos(Angle+4.0/5.0*PI-1.0/15.0*PI);
	P[3][4]=2.0/5.0*cos(Angle+2.0/5.0*PI-1.0/15.0*PI);

	P[4][0]=-2.0/5.0*sin(Angle-1.0/15.0*PI);
	P[4][1]=-2.0/5.0*sin(Angle-2.0/5.0*PI-1.0/15.0*PI);
	P[4][2]=-2.0/5.0*sin(Angle-4.0/5.0*PI-1.0/15.0*PI);
	P[4][3]=-2.0/5.0*sin(Angle+4.0/5.0*PI-1.0/15.0*PI);
	P[4][4]=-2.0/5.0*sin(Angle+2.0/5.0*PI-1.0/15.0*PI);

	P[5][0]=1.0/5.0;
	P[5][1]=1.0/5.0;
	P[5][2]=1.0/5.0;
	P[5][3]=1.0/5.0;
	P[5][4]=1.0/5.0;

	P[6][0]=2.0/5.0*cos(Angle-2.0/15.0*PI);
	P[6][1]=2.0/5.0*cos(Angle-2.0/5.0*PI-2.0/15.0*PI);
	P[6][2]=2.0/5.0*cos(Angle-4.0/5.0*PI-2.0/15.0*PI);
	P[6][3]=2.0/5.0*cos(Angle+4.0/5.0*PI-2.0/15.0*PI);
	P[6][4]=2.0/5.0*cos(Angle+2.0/5.0*PI-2.0/15.0*PI);

	P[7][0]=-2.0/5.0*sin(Angle-2.0/15.0*PI);
	P[7][1]=-2.0/5.0*sin(Angle-2.0/5.0*PI-2.0/15.0*PI);
	P[7][2]=-2.0/5.0*sin(Angle-4.0/5.0*PI-2.0/15.0*PI);
	P[7][3]=-2.0/5.0*sin(Angle+4.0/5.0*PI-2.0/15.0*PI);
	P[7][4]=-2.0/5.0*sin(Angle+2.0/5.0*PI-2.0/15.0*PI);

	P[8][0]=1.0/5.0;
	P[8][1]=1.0/5.0;
	P[8][2]=1.0/5.0;
	P[8][3]=1.0/5.0;
	P[8][4]=1.0/5.0;
}
void Motor15::invParkTransMatrix(double INRP[15][3],double Angle)
{
	INRP[0][0]=cos(Angle);
	INRP[1][0]=cos(Angle-2.0/5.0*PI);
	INRP[2][0]=cos(Angle-4.0/5.0*PI);
	INRP[3][0]=cos(Angle+4.0/5.0*PI);
	INRP[4][0]=cos(Angle+2.0/5.0*PI);

	INRP[0][1]=-sin(Angle);
	INRP[1][1]=-sin(Angle-2.0/5.0*PI);
	INRP[2][1]=-sin(Angle-4.0/5.0*PI);
	INRP[3][1]=-sin(Angle+4.0/5.0*PI);
	INRP[4][1]=-sin(Angle+2.0/5.0*PI);

	INRP[0][2]=1.0;
	INRP[1][2]=1.0;
	INRP[2][2]=1.0;
	INRP[3][2]=1.0;
	INRP[4][2]=1.0;

	INRP[5][0]=cos(Angle-1.0/15.0*PI);
	INRP[6][0]=cos(Angle-2.0/5.0*PI-1.0/15.0*PI);
	INRP[7][0]=cos(Angle-4.0/5.0*PI-1.0/15.0*PI);
	INRP[8][0]=cos(Angle+4.0/5.0*PI-1.0/15.0*PI);
	INRP[9][0]=cos(Angle+2.0/5.0*PI-1.0/15.0*PI);

	INRP[5][1]=-sin(Angle-1.0/15.0*PI);
	INRP[6][1]=-sin(Angle-2.0/5.0*PI-1.0/15.0*PI);
	INRP[7][1]=-sin(Angle-4.0/5.0*PI-1.0/15.0*PI);
	INRP[8][1]=-sin(Angle+4.0/5.0*PI-1.0/15.0*PI);
	INRP[9][1]=-sin(Angle+2.0/5.0*PI-1.0/15.0*PI);

	INRP[5][2]=1.0;
	INRP[6][2]=1.0;
	INRP[7][2]=1.0;
	INRP[8][2]=1.0;
	INRP[9][2]=1.0;

	INRP[10][0]=cos(Angle-2.0/15.0*PI);
	INRP[11][0]=cos(Angle-2.0/5.0*PI-2.0/15.0*PI);
	INRP[12][0]=cos(Angle-4.0/5.0*PI-2.0/15.0*PI);
	INRP[13][0]=cos(Angle+4.0/5.0*PI-2.0/15.0*PI);
	INRP[14][0]=cos(Angle+2.0/5.0*PI-2.0/15.0*PI);

	INRP[10][1]=-sin(Angle-2.0/15.0*PI);
	INRP[11][1]=-sin(Angle-2.0/5.0*PI-2.0/15.0*PI);
	INRP[12][1]=-sin(Angle-4.0/5.0*PI-2.0/15.0*PI);
	INRP[13][1]=-sin(Angle+4.0/5.0*PI-2.0/15.0*PI);
	INRP[14][1]=-sin(Angle+2.0/5.0*PI-2.0/15.0*PI);

	INRP[10][2]=1.0;
	INRP[11][2]=1.0;
	INRP[12][2]=1.0;
	INRP[13][2]=1.0;
	INRP[14][2]=1.0;
}

//�ж��Ƿ��нڵ�ӵ�
bool Motor15::isSeries(int from,int to)
{
	int series = 1;

	if(to==0) {to=from;series=0;}
	if(from==0) {from=to;series=0;}
	if(to==0) cerr<<"Both Nodes Zero!!"<<endl;

	return series;
}