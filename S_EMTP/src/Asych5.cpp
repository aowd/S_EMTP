#include "Asych5.h"

#include <cmath>
#include <iostream>
using namespace std;

Asych5::Asych5(int firstNode,int lastNode)
{
	//�Զ���Ԫ��
	isUserDef = 1;
	need_NEC=1;
	type=555;
	nPort=10;
	nodeNumber=new int[nPort];
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}

	//dq0����ϵ�µĵ������
	Rp = 0.0086917;
	Rr1 = 0.0047;
	Rrt = 0.161;
	//dq0����ϵ�µĵ翹����
	Xl1 = 0.202625;
	Xlr1 = 0.077525;
	Xm1 = 1.7819;
	Xlt1 = 0.01;
	Xlrt = 0.1492;
	Xmt = 8.117;
	Xs0 = 0.15045;
	//�м����
	Xs1 = Xl1 + Xm1; 
	Xr1 = Xlr1 + Xm1;
	Xst = Xlt1 + Xmt;
	Xrt = Xlrt + Xmt;

	//���Ӳ��ѹ�������ͽ�Ƶ�ʵĻ�ֵ
	Vbase = 2400;//V,peak value
	Ibase = 2052;//A,peak value
	wbase = 2*PI*20;//fbase=20Hz;

	//��е�����
	Je = 25000;//ת������


	//�ɿ˱任�����������
	P.ResizeTo(1,3,1,5);
	invP.ResizeTo(1,5,1,3);

	/*     ���λ��ַ���Ҫʹ�õ�ϵ������      */
	coff = 39.0/41.0;//alpha
	T=1.0/(deltaT*wbase);
	//
	A.ResizeTo(1,5,1,5);//(5,5)
	B.ResizeTo(1,5,1,5);//(5,5)
	AA.ResizeTo(1,5,1,5);//(5,5)
	BB.ResizeTo(1,5,1,5);//(5,5)
	//
	Rss.ResizeTo(1,3,1,3);//AA(1:3,1:3)
	Rsr.ResizeTo(1,3,1,2);//AA(1:3,4:5)
	Rrs.ResizeTo(1,2,1,3);//AA(4:5,1:3)
	Rrr.ResizeTo(1,2,1,2);//AA(4:5,4:5)
	//
	Bss.ResizeTo(1,3,1,3);//BB(1:3,1:3)
	Bsr.ResizeTo(1,3,1,2);//BB(1:3,4:5)
	Brs.ResizeTo(1,2,1,3);//BB(4:5,1:3)
	Brr.ResizeTo(1,2,1,2);//BB(4:5,4:5)

	//����ϵ������ʹ��ת����abcde����ϵʱ�ܽ���
	Rs.ResizeTo(1,3,1,3);//
	Rs1.ResizeTo(1,3,1,3);//
	Rs2.ResizeTo(1,3,1,3);//Rs2=Rs-Rs1;

	//
	Gs1.ResizeTo(1,3,1,3);//
	Gs2.ResizeTo(1,3,1,3);//
	Gs3.ResizeTo(1,3,1,3);//
	Gs4.ResizeTo(1,3,1,2);//

	//
	Gr1.ResizeTo(1,2,1,2);//
	Gr2.ResizeTo(1,2,1,3);//
	Gr3.ResizeTo(1,2,1,3);//
	Gr4.ResizeTo(1,2,1,2);//

	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=0;i<10;i++)
	{
		Vnode[i] = 0;//�ڵ��ѹ
	}
	for (int i=0;i<5;i++)
	{
		Vabcde[i] = 0;//֧·��ѹ
		Iabcde[i] = 0;//֧·����
		Vabcde_1[i] = 0;
		Iabcde_1[i] = 0;
		Vabcde_2[i] = 0;
		Iabcde_2[i] = 0;
	}

	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=0;i<3;i++)
	{
		Vdqs[i] = 0;
		Idqs[i] = 0;
		Idqs_his[i] = 0;//Ԥ����
		Idqs_his2[i] = 0;//Ԥ����
		Idqs_forcast[i] = 0;//Ԥ��ֵ
	}

	//dq0����ϵ�¶��Ӳ�Ĵ���
	psi.ResizeTo(1,5);
	Isr.ResizeTo(1,5);//�������ʱ�õĵ���

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
	for (int i=0;i<5;i++)
	{
		nortonEquivalentCurrent[i] = 0;//abcde�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_1[i] = 0;
		nortonEquivalentCurrent_2[i] = 0;
	}
	for (int i=0;i<3;i++)
	{
		nortonEquivalentCurrent_dq0[i] = 0;//dq0�����µ�ŵ�ٵ�Ч����������ֵ
	}

	//����Gs2,Gs3,Gs4
	calculateCoefficientMatrix();
}

//��������
Asych5::~Asych5()
{
}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/
//��ʼ��֧·��ѹ����
void Asych5::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	for (int k=0;k<5;k++)
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
	for (int i=0;i<3;i++)
	{
		Vdqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Vdqs[i] += P(i+1,j+1)*Vabcde[j];
		}
	}

	//����dq0����ϵ�µĶ��ӵ���
	for (int i=0;i<3;i++)
	{
		Idqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Idqs[i] += P(i+1,j+1)*Iabcde[j];
		}
	}

	for (int k=0;k<3;k++)
	{
		Idqs_his2[k] = Idqs[k];
		Idqs_his[k] = Idqs[k];
	}

	//Ԥ����һʱ�̶��ӵ���
	for (int i=0;i<3;i++)
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
		for (int j=0;j<3;j++)
		{
			Ir[i] = Ir[i] + Gr2(i+1,j+1)*Idqs[j] + Gr3(i+1,j+1)*Idqs_his[j];
		}
		Ir[i] = Ir[i] + Gr4(i+1,1)*Ir_his[0] + Gr4(i+1,2)*Ir_his[1];
	}

	//����Inew
	TVectorD tem1,Idqs_forcast_V,Idqs_V,Vdqs_V,Ir_V,Inew;
	tem1.ResizeTo(1,3);
	Idqs_forcast_V.ResizeTo(1,3);
	Idqs_V.ResizeTo(1,3);
	Vdqs_V.ResizeTo(1,3);
	Ir_V.ResizeTo(1,2);
	Inew.ResizeTo(1,3);

	for (int i=1;i<=3;i++)
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
	for (int i=0;i<3;i++)
	{
		nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	}

	Angle = wr*deltaT + theta;
	invParkTransMatrix(invP,Angle);

	for (int i=0;i<5;i++)
	{
		nortonEquivalentCurrent[i] = 0;
		for (int j=0;j<3;j++)
		{
			nortonEquivalentCurrent[i] += invP(i+1,j+1)*nortonEquivalentCurrent_dq0[j];
		}
	}

	//���¼���Idqs��Idqs_his��Idqs_his2
	for (int k=0;k<3;k++)
	{
		Idqs[k] = Inew(k+1);
		Idqs_his[k] = Idqs[k];
		Idqs_his2[k] = Idqs[k];
	}

	//����е����
	//����Isr
	for (int i=1;i<=3;i++)
	{
		Isr(i) = Idqs[i-1];
	}
	for (int i=1;i<=2;i++)
	{
		Isr(i+3) = Ir[i-1];
	}

	//�������psi
	psi = B*Isr;

	for (int i=1;i<=5;i++)
	{
		psi(i) /= wbase;
	}

	//������ת��Te
	Te = (psi(1)*Idqs[1]-psi(2)*Idqs[0])*15.0/1e6;

	//����ת�ӽ�Ƶ��wr
	wr_his = wr;
	wr = wr_his + (Te-Tm)*deltaT/Je;

	//����ת�ӽ�theta
	theta = theta+(wr-wr0)*deltaT;

	//��PSCAD���ݳ�ʼ�����ӵ�ѹ����
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<5;i++)
	{
		Iabcde[i] = initialCurrentArray[ptr+i];
	}
	ptr+=5;
}

//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
void Asych5::readNodeVoltage(TVectorD& nodeVoltageArray)
{//�ӽڵ��ѹ�����ж��綯���˽ڵ��ѹ
	for (int i=0;i<10;i++)
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
void Asych5::calculateBranchVoltage()
{
	for (int i=0;i<5;i++)
	{
		Vabcde_1[i] = Vabcde[i];
		Vabcde[i] = Vnode[2*i] - Vnode[2*i+1];
	}
}

//����֧·����
void Asych5::calculateBranchCurrent()
{
	for (int i=0;i<5;i++)
	{
		Iabcde_1[i] = Iabcde[i];
		Iabcde[i] = nortonEquivalentConductance*Vabcde[i] + nortonEquivalentCurrent[i];
	}

	//cout<<endl;
	//for (int i=0;i<5;i++)
	//{
	//	cout<<Vabcde[i]<<endl;
	//}
}

//����֧·��ŵ�ٵ�Ч��·�еĵ�����
void Asych5::calculateNortonEquivalentCurrent(double time)
{
	calculateCoefficientMatrix();
	calculateDq0Results(time);

	TVectorD tem1,Idqs_forcast_V,Idqs_V,Vdqs_V,Ir_V,Inew;
	tem1.ResizeTo(1,3);
	Idqs_forcast_V.ResizeTo(1,3);
	Idqs_V.ResizeTo(1,3);
	Vdqs_V.ResizeTo(1,3);
	Ir_V.ResizeTo(1,2);
	Inew.ResizeTo(1,3);

	for (int i=1;i<=3;i++)
	{
		Idqs_forcast_V(i) = Idqs_forcast[i-1];
		Idqs_V(i) = Idqs[i-1];
		Vdqs_V(i) = Vdqs[i-1];
	}
	Ir_V(1) = Ir[0];
	Ir_V(2) = Ir[1];

	tem1 = Gs2*Idqs_forcast_V + Gs3*Idqs_V + Gs4*Ir_V + coff*Vdqs_V;
	Inew = Gs1*tem1;

	for (int i=0;i<3;i++)
	{
		nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	}

	//cout<<endl<<"Idqs_forcast_V:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Idqs_forcast_V(i)<<endl;
	//}

	//cout<<endl<<"Idqs_V:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Idqs_V(i)<<endl;
	//}

	//cout<<endl<<"Ir_V:"<<endl;
	//for (int i=1;i<=2;i++)
	//{
	//	cout<<Ir_V(i)<<endl;
	//}

	//cout<<endl<<"Vdqs_V:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Vdqs_V(i)<<endl;
	//}

	//cout<<endl<<"Inew:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Inew(i)<<endl;
	//}

	double Angle;
	Angle = wr0*time + wr*deltaT + theta;
	invParkTransMatrix(invP,Angle);

	for (int k=0;k<5;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	for (int i=0;i<5;i++)
	{
		nortonEquivalentCurrent[i] = 0;
		for (int j=0;j<3;j++)
		{
			nortonEquivalentCurrent[i] += invP(i+1,j+1)*nortonEquivalentCurrent_dq0[j];
		}
	}

	//cout<<endl<<"nortonEquivalentCurrent:"<<endl;
	//for (int i=0;i<5;i++)
	//{
	//	cout<<nortonEquivalentCurrent[i]<<endl;
	//}

	calculateOmega();
}

//�γɽڵ�ŵ�ٵ�Ч��������
void Asych5::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{
	int k,m;
	for (int i=0;i<5;i++)
	{
		k = nodeNumber[2*i];
		m = nodeNumber[2*i+1];

		nodeNortonEquivalentCurrentArray(k) -= nortonEquivalentCurrent[i];
		nodeNortonEquivalentCurrentArray(m) += nortonEquivalentCurrent[i];
	}
}
void Asych5::formConductanceMatrix(TMatrixD &conductanceMatrix)
{
	int k,m;
	for (int i=0;i<5;i++)
	{
		k = nodeNumber[2*i];
		m = nodeNumber[2*i+1];

		conductanceMatrix(k,k) += nortonEquivalentConductance;
		conductanceMatrix(m,m) += nortonEquivalentConductance;
		conductanceMatrix(k,m) -= nortonEquivalentConductance;
		conductanceMatrix(m,k) -= nortonEquivalentConductance;
	}
}

//����֧·����
void Asych5::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{
	for (int i=0;i<5;i++)
	{
		branchCurrentMatrix(counter,ptr+i) = Iabcde[i];
	}
	ptr+=5;
}

//����֧·����
void Asych5::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int i=0;i<5;i++)
	{
		branchCurrentMatrix_1[counter-1][ptr+i] = Iabcde[i];
	}
	ptr+=5;
}

//������ֵ����֧·�ĵ�ѹ�������в�ֵ
void Asych5::interpolate(double ratio)
{
	//��ѹ����������ֵ
	for (int k=0;k<5;k++)
	{
		Vabcde[k] = (1-ratio)*Vabcde_1[k] + ratio*Vabcde[k];
		Iabcde[k] = (1-ratio)*Iabcde_1[k] + ratio*Iabcde[k];
	}

	//ŵ�ٵ�ֵ������ֵ
	double temp1,temp2;
	for (int k=0;k<5;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void Asych5::updateResult(int updateTypeNum)
{
	double V_temp[5],I_temp[5];//��ʱ��������������ʱʹ��

	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for (int i=0;i<5;i++)
		{
			Vabcde_2[i] = Vabcde_1[i];
			Iabcde_2[i] = Iabcde_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://��_2��������ֵ����_1������
		for (int i=0;i<5;i++)
		{
			Vabcde_1[i] = Vabcde_2[i];
			Iabcde_1[i] = Iabcde_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://ŵ�ٵ�ֵ�����洢����_1��������ֵ����_2������
		for (int i=0;i<5;i++)
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
void Asych5::calculateCoefficientMatrix()
{
	double w;
	w=wr/wbase;

	//����ϵ������A��B
	A(1,1)=Rp;
	A(1,2)=-w*Xs1;
	A(1,5)=-w*Xm1;

	A(2,1)=w*Xs1;
	A(2,2)=Rp;
	A(2,4)=w*Xm1;

	A(3,3)=Rp;

	A(4,4)=Rr1;
	A(5,5)=Rr1;

	B(1,1)=Xs1;
	B(1,4)=Xm1;

	B(2,2)=Xs1;
	B(2,5)=Xm1;

	B(3,3)=Xs0;

	B(4,1)=Xm1;
	B(4,4)=Xr1;

	B(5,2)=Xm1;
	B(5,5)=Xr1;

	//����ϵ������AA��BB
	AA=A+(1.0+coff)*T*B;
	BB=coff*A-(1.0+coff)*T*B;

	//����ϵ������Rss...Brr
	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rss(i,j) = AA(i,j);
			Bss(i,j) = BB(i,j);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rsr(i,j) = AA(i,j+3);
			Bsr(i,j) = BB(i,j+3);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rrs(i,j) = AA(i+3,j);
			Brs(i,j) = BB(i+3,j);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rrr(i,j) = AA(i+3,j+3);
			Brr(i,j) = BB(i+3,j+3);
		}
	}

	//����ϵ������Gr1...Gr4
	Gr1(1,1)=-1.0/Rrr(1,1);
	Gr1(2,2)=Gr1(1,1);

	Gr2 = Gr1*Rrs;
	Gr3 = Gr1*Brs;
	Gr4 = Gr1*Brr;

	//����ϵ������Reqs...inv_Rr,Gs1...Gs4
	Rs = Rss+Rsr*Gr2;
	Tp = Rs(1,1);

	Rs1(1,1)=Tp;
	Rs1(2,2)=Tp;
	Rs1(3,3)=Tp;

	Rs2 = Rs-Rs1;
 
	nortonEquivalentConductance = 1/Tp;
	//cout<<nortonEquivalentConductance<<endl;

	Gs1(1,1)=1.0/Rs(1,1);
	Gs1(2,2)=Gs1(1,1);
	Gs1(3,3)=Gs1(1,1);

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Gs2(i,j) = -1*Rs2(i,j);
		}
	}

	Gs3 = Bss+Rsr*Gr3;
	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Gs3(i,j) = -1*Gs3(i,j);
		}
	}

	Gs4 = Bsr+Rsr*Gr4;
	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Gs4(i,j) = -1*Gs4(i,j);
		}
	}
}

//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
void Asych5::calculateDq0Results(double time)
{
	//����Park�任����
	double Angle;
	Angle = wr0*time+theta;
	parkTransMatrix(P,Angle);

	//����dq0����ϵ�µĶ��ӵ�ѹ
	for (int i=0;i<3;i++)
	{
		Vdqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Vdqs[i] += P(i+1,j+1)*Vabcde[j];
		}
	}

	//����dq0����ϵ�µĶ��ӵ���
	for (int i=0;i<3;i++)
	{
		Idqs_his2[i] = Idqs_his[i];
		Idqs_his[i] = Idqs[i];
		Idqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Idqs[i] += P(i+1,j+1)*Iabcde[j];
		}
	}

	//Ԥ����һʱ�̶��ӵ���
	for (int i=0;i<3;i++)
	{
		Idqs_forcast[i] = (202.0/120.0)*Idqs[i]-(44.0/120.0)*Idqs_his[i]-(38.0/120.0)*Idqs_his2[i];
	}

	//����ת�ӵ���
	Ir_his[0] = Ir[0];
	Ir_his[1] = Ir[1];

	for (int i=0;i<2;i++)
	{
		Ir[i] = 0;
		for (int j=0;j<3;j++)
		{
			Ir[i] = Ir[i] + Gr2(i+1,j+1)*Idqs[j] + Gr3(i+1,j+1)*Idqs_his[j];
		}
		Ir[i] = Ir[i] + Gr4(i+1,1)*Ir_his[0] + Gr4(i+1,2)*Ir_his[1];
	}

	//for (int i=0;i<2;i++)
	//{
	//	cout<<Ir[i]<<endl;
	//}
}

//����е���̣�����ת��ת��
void Asych5::calculateOmega()
{
	////����Isr
	for (int i=1;i<=3;i++)
	{
		Isr(i) = Idqs[i-1];
	}
	for (int i=1;i<=2;i++)
	{
		Isr(i+3) = Ir[i-1];
	}

	//�������psi
	psi = B*Isr;

	for (int i=1;i<=5;i++)
	{
		psi(i) /= wbase;
	}

	//������ת��Te
	Te = (psi(1)*Idqs[1]-psi(2)*Idqs[0])*15.0/1e6;

	//����ת�ӽ�Ƶ��wr
	wr_his = wr;
	wr = wr_his + (Te-Tm)*deltaT/Je;

	//����ת�ӽ�theta
	theta = theta+(wr-wr0)*deltaT;

	//cout<<wr<<endl;
	//cout<<theta<<endl;
}

//�����ɿ˱任����
void Asych5::parkTransMatrix(TMatrixD &P,double Angle)
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
}

//�����ɿ˷��任����
void Asych5::invParkTransMatrix(TMatrixD &INRP,double Angle)
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
}


//�ж��Ƿ��нڵ�ӵ�
bool Asych5::isSeries(int from,int to)
{
	int series = 1;

	if(to==0) {to=from;series=0;}
	if(from==0) {from=to;series=0;}
	if(to==0) cerr<<"Both Nodes Zero!!"<<endl;

	return series;
}