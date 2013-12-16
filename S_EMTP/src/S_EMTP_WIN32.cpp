#ifdef WIN32
#include "EMTP.h"
void S_EMTP_WIN32()
{
	EMTP emtp;
	emtp.specifySystem();//����ϵͳ���뺯����ʹ��Ĭ�ϲ��������沽��50us��
	emtp.specifyMsrSystem();//����ϵͳ���뺯��
	emtp.specifyCtrlSystem();//����ϵͳ���뺯��

	// ��ʼ��
	emtp.counter = 1;
	emtp.counter2 = 0;
	emtp.cnt_p = 0;
	emtp.curTime = emtp.startTime;
	emtp.timeArray(1) = emtp.startTime;
	emtp.initializeSystem();

	//===================��ʱ��ʼ===================
	int nstep = 0; 
	cout << "Run begins!" << endl;
	emtp.my_clock(emtp.tStartAll);
	
	double* tau_new = new double[3*emtp.nPWMConverter];
	int correctOrNot;
	// ��ѭ��
	while ( emtp.cnt_p < emtp.Np )
	{
		//cout<<emtp.curTime<<endl;
	    emtp.solveG(tau_new);
		correctOrNot = emtp.check_tau(tau_new,0.005);
		if (correctOrNot == 1)
		{
			emtp.correct_tau(tau_new);
			emtp.restore_state();
		} 
		else
		{
			emtp.cnt_p++;
			if (emtp.cnt_p < emtp.Np)
			{
				emtp.predict_tau();
				emtp.store_state();
			}
		}
	}
	emtp.my_clock(emtp.tEndAll);
	cout << "Run ends!" << endl;

	// ��������
	//emtp.saveResult(); //ϵͳ�������������������2ά����ĸ�ʽ�洢���ı��ļ�result.dat��
	delete emtp.lu_conductanceMatrix;//ÿһ��case����new������lu����Ҫɾ����
	cout<<"�ܺ�ʱΪ\t"<<emtp.cal_timeCost(emtp.tStartAll,emtp.tEndAll)/1e6<<"s"<<endl;
	cout<<"ÿ����ƽ����ʱΪ\t"<<emtp.cal_timeCost(emtp.tStartAll,emtp.tEndAll)/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	cout << "ƽ����������Ϊ\t" <<emtp.nSolveG/emtp.Np<< endl;
	//cout<<"����Ԫ��ŵ�ٵ�ֵ����\t"<<emtp.t1/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	//cout<<"�γɽڵ�ע���������\t"<<emtp.t2/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	//cout<<"���ڵ��ѹ����\t"<<emtp.t3/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	//cout<<"����Ԫ���ڲ�����\t"<<emtp.t4/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	cout<<"����ϵͳ���\t"<<emtp.t5/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	cout<<"����ϵͳ���\t"<<emtp.t6/emtp.nSolveG/emtp.Ns<<"us"<<endl;

	//		EMTP emtp;
//
//		emtp.specifySystem();//����ϵͳ���뺯����ʹ��Ĭ�ϲ��������沽��50us��
//		emtp.specifyMsrSystem();//����ϵͳ���뺯��
//		emtp.specifyCtrlSystem();//����ϵͳ���뺯��
//
//		// ������ѭ��
//		emtp.counter = 1;//������
//		emtp.counter2 = 0;
//		emtp.curTime = emtp.startTime;//��ǰʱ��
//		emtp.timeArray(1) = emtp.startTime;
//
//		emtp.initializeSystem();//��ʼ��ϵͳ��0--0������1--PSCAD��ֵ����
//
//		// Update current time and counter
//		emtp.curTime= emtp.curTime + emtp.deltaT;
//		emtp.counter =emtp.counter + 1;
//		emtp.counter2 += 1;
//		emtp.timeArray(emtp.counter) = emtp.curTime;
//		
//		emtp.transferControlVariables();//�������ź�ֵ����������ϵͳ
//		emtp.switchTreatment();// Check switch state
//		emtp.calculateBranchNortonEquivalentCurrent(emtp.curTime);//����֧·ŵ�ٵ�Ч����
//		emtp.formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������
//
//#ifndef PREDICT
//		emtp.formConductanceMatrix();
//#endif
//
////===================��ʱ��ʼ===================
//		int nstep = 0; 
//		cout << "Run begins!" << endl;
//		clock_t t1 = clock();	
//		while(emtp.curTime<emtp.finishTime-emtp.deltaT/2)
//		{
//			nstep++;
//			
//			emtp.solveNodeVoltageEquation();//���ڵ��ѹ����
//			emtp.saveNodeVoltage();//����ڵ��ѹ����,�����������ȫ�ֵ�ѹ�����У����֧·�޹ء�
//			emtp.calculateBranchVoltage();//����֧·��ѹ������
//			emtp.calculateBranchCurrent();//����֧·����
//			emtp.saveBranchCurrent();//����֧·����			
//			emtp.transferMeasurands(); // ������ϵͳ����ֵ���䵽����ϵͳ
//			emtp.solveCtrlSystem();
//
//			/////////////////////////////////////////////////////////////////////////////
//			// У�����ض���ʱ��
//			if (emtp.nPWMConverter != 0) {
//#ifdef GENERAL
//				emtp.correctSwitchingInstants(0.005);
//#else
//				emtp.correctSwitchingInstants(0.005);
//#endif
//			}
//			/////////////////////////////////////////////////////////////////////////////
//
//			emtp.curTime= emtp.curTime + emtp.deltaT;
//			emtp.counter =emtp.counter + 1;
//			emtp.counter2 += 1;
//			emtp.timeArray(emtp.counter) = emtp.curTime;
//	
//			emtp.transferControlVariables(); // �������ź�ֵ����������ϵͳ		
//			emtp.switchTreatment();// Check switch state
//			emtp.calculateBranchNortonEquivalentCurrent(emtp.curTime);//����֧·ŵ�ٵ�Ч����
//			emtp.formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������
//
//#ifndef PREDICT
//		    emtp.formConductanceMatrix();
//#endif
//
//		}
//
//		emtp.solveNodeVoltageEquation();//���ڵ��ѹ����
//		emtp.saveNodeVoltage();//����ڵ��ѹ����,�����������ȫ�ֵ�ѹ�����У����֧·�޹ء�
//		emtp.calculateBranchVoltage();//����֧·��ѹ������
//		emtp.calculateBranchCurrent();//����֧·����
//		emtp.saveBranchCurrent();//����֧·����
//		emtp.transferMeasurands(); // ������ϵͳ����ֵ���䵽����ϵͳ
//		emtp.solveCtrlSystem();
//		emtp.switchTreatment();// Check switch state
//      
//		clock_t t2 = clock();	
//		cout << "Run ends!" << endl;
////===================��ʱ����===================			
//		emtp.saveResult(); //ϵͳ�������������������2ά����ĸ�ʽ�洢���ı��ļ�result.dat��
////		emtp.checkResult();
//		delete emtp.lu_conductanceMatrix;//ÿһ��case����new������lu����Ҫɾ����
//        
//		nstep = ( emtp.finishTime - emtp.startTime ) / emtp.deltaT;
//		cout<<"clock time = "<<(double)(t2-t1)/nstep*1000<<"us"<<endl;
//		cout<<"�ܺ�ʱΪ"<<(double)(t2-t1)/1000<<"s"<<endl;
}
#endif
