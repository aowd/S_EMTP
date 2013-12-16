#ifdef WIN32
#include "EMTP.h"
void S_EMTP_WIN32()
{
	EMTP emtp;
	emtp.specifySystem();//电气系统输入函数：使用默认参数，仿真步长50us。
	emtp.specifyMsrSystem();//测量系统输入函数
	emtp.specifyCtrlSystem();//控制系统输入函数

	// 初始化
	emtp.counter = 1;
	emtp.counter2 = 0;
	emtp.cnt_p = 0;
	emtp.curTime = emtp.startTime;
	emtp.timeArray(1) = emtp.startTime;
	emtp.initializeSystem();

	//===================计时开始===================
	int nstep = 0; 
	cout << "Run begins!" << endl;
	emtp.my_clock(emtp.tStartAll);
	
	double* tau_new = new double[3*emtp.nPWMConverter];
	int correctOrNot;
	// 主循环
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

	// 后续处理
	//emtp.saveResult(); //系统输出函数：将计算结果按2维数组的格式存储在文本文件result.dat中
	delete emtp.lu_conductanceMatrix;//每一个case里面new出来的lu对象都要删除掉
	cout<<"总耗时为\t"<<emtp.cal_timeCost(emtp.tStartAll,emtp.tEndAll)/1e6<<"s"<<endl;
	cout<<"每步长平均耗时为\t"<<emtp.cal_timeCost(emtp.tStartAll,emtp.tEndAll)/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	cout << "平均迭代次数为\t" <<emtp.nSolveG/emtp.Np<< endl;
	//cout<<"计算元件诺顿等值电流\t"<<emtp.t1/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	//cout<<"形成节点注入电流向量\t"<<emtp.t2/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	//cout<<"求解节点电压方程\t"<<emtp.t3/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	//cout<<"计算元件内部变量\t"<<emtp.t4/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	cout<<"电气系统求解\t"<<emtp.t5/emtp.nSolveG/emtp.Ns<<"us"<<endl;
	cout<<"控制系统求解\t"<<emtp.t6/emtp.nSolveG/emtp.Ns<<"us"<<endl;

	//		EMTP emtp;
//
//		emtp.specifySystem();//电气系统输入函数：使用默认参数，仿真步长50us。
//		emtp.specifyMsrSystem();//测量系统输入函数
//		emtp.specifyCtrlSystem();//控制系统输入函数
//
//		// 进入主循环
//		emtp.counter = 1;//计数器
//		emtp.counter2 = 0;
//		emtp.curTime = emtp.startTime;//当前时间
//		emtp.timeArray(1) = emtp.startTime;
//
//		emtp.initializeSystem();//初始化系统，0--0启动；1--PSCAD初值启动
//
//		// Update current time and counter
//		emtp.curTime= emtp.curTime + emtp.deltaT;
//		emtp.counter =emtp.counter + 1;
//		emtp.counter2 += 1;
//		emtp.timeArray(emtp.counter) = emtp.curTime;
//		
//		emtp.transferControlVariables();//将控制信号值传输至电气系统
//		emtp.switchTreatment();// Check switch state
//		emtp.calculateBranchNortonEquivalentCurrent(emtp.curTime);//计算支路诺顿等效电流
//		emtp.formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量
//
//#ifndef PREDICT
//		emtp.formConductanceMatrix();
//#endif
//
////===================计时开始===================
//		int nstep = 0; 
//		cout << "Run begins!" << endl;
//		clock_t t1 = clock();	
//		while(emtp.curTime<emtp.finishTime-emtp.deltaT/2)
//		{
//			nstep++;
//			
//			emtp.solveNodeVoltageEquation();//求解节点电压方程
//			emtp.saveNodeVoltage();//保存节点电压向量,将求解结果存在全局电压矩阵中，与各支路无关。
//			emtp.calculateBranchVoltage();//计算支路电压并保存
//			emtp.calculateBranchCurrent();//计算支路电流
//			emtp.saveBranchCurrent();//保存支路电流			
//			emtp.transferMeasurands(); // 将电气系统测量值传输到控制系统
//			emtp.solveCtrlSystem();
//
//			/////////////////////////////////////////////////////////////////////////////
//			// 校正开关动作时刻
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
//			emtp.transferControlVariables(); // 将控制信号值传输至电气系统		
//			emtp.switchTreatment();// Check switch state
//			emtp.calculateBranchNortonEquivalentCurrent(emtp.curTime);//计算支路诺顿等效电流
//			emtp.formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量
//
//#ifndef PREDICT
//		    emtp.formConductanceMatrix();
//#endif
//
//		}
//
//		emtp.solveNodeVoltageEquation();//求解节点电压方程
//		emtp.saveNodeVoltage();//保存节点电压向量,将求解结果存在全局电压矩阵中，与各支路无关。
//		emtp.calculateBranchVoltage();//计算支路电压并保存
//		emtp.calculateBranchCurrent();//计算支路电流
//		emtp.saveBranchCurrent();//保存支路电流
//		emtp.transferMeasurands(); // 将电气系统测量值传输到控制系统
//		emtp.solveCtrlSystem();
//		emtp.switchTreatment();// Check switch state
//      
//		clock_t t2 = clock();	
//		cout << "Run ends!" << endl;
////===================计时结束===================			
//		emtp.saveResult(); //系统输出函数：将计算结果按2维数组的格式存储在文本文件result.dat中
////		emtp.checkResult();
//		delete emtp.lu_conductanceMatrix;//每一个case里面new出来的lu对象都要删除掉
//        
//		nstep = ( emtp.finishTime - emtp.startTime ) / emtp.deltaT;
//		cout<<"clock time = "<<(double)(t2-t1)/nstep*1000<<"us"<<endl;
//		cout<<"总耗时为"<<(double)(t2-t1)/1000<<"s"<<endl;
}
#endif
