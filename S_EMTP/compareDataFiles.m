clc
clear all;
close all;
format long
stg=strcat('================',datestr(now),'================');
disp(stg);
%% ������Ϣ
plotFlag = 1;%��ͼ��־λ
iCase=435;%�������
switch iCase
    case 211
        nNode = 15;%�ڵ���
        loadID = [1 4 1-nNode 4-nNode];
        nFile = 4;
    case 213
        nNode = 12;%�ڵ���
        loadID = [1 4 1-nNode 4-nNode];
        nFile = 3;
     case 214
        nNode = 12;%�ڵ���
        loadID = [1 4 1-nNode 4-nNode];
        nFile = 3;
    case 215
        nNode = 6;%�ڵ���
        loadID = [1 4 1-nNode 4-nNode];
        nFile = 2;
    case 272
        nNode = 5;%�ڵ�������������������
        loadID = [1 4 1-nNode 4-nNode];
        nFile = 2;
    case 2721
        nNode = 12;%�ڵ���
        loadID = [1 4 1-nNode 4-nNode];
        nFile = 4;
    case 291
        nNode = 11;%�ڵ���
        loadID = [1-nNode 4-nNode 1 7];
        nFile = 3;
    case 292
        nNode = 11;%�ڵ���
        loadID = [1-nNode 4-nNode 1 7];
        nFile = 3;
    case 294
        nNode = 10;%�ڵ���
        loadID = [1-nNode 4-nNode 1 7];
        nFile = 3;
    case 295
        nNode = 9;%�ڵ���
        loadID = [1-nNode 4-nNode 1 7];
        nFile = 3;
    case 407
        nNode = 17;%�ڵ���
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 5;
     case 408
        nNode = 17;%�ڵ���
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 5;
     case 4082
        nNode = 17*2;%�ڵ���
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 9;
     case 409
        nNode = 136;%�ڵ���
        loadID = [39-nNode 32-nNode 60 53];
        nFile = 35;
     case 411
        nNode = 10;%�ڵ���
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 3;
     case 412
        nNode = 66;%�ڵ���
        loadID = [1-nNode 32-nNode 1 53];
        nFile = 17;
     case 414
        nNode = 21;%�ڵ���
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 6;
     case 415
        nNode = 23;%�ڵ���
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 6;
     case 417
        nNode = 73;%�ڵ����������������
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 21;
     case 4171
        nNode = 17;%�ڵ�������糡ϵ������
        loadID = [1-nNode 4-nNode 1 4];
        nFile = 5;
     case 420
        nNode = 5;%�ڵ������������ģ��
        loadID = [1-nNode 2-nNode 3-nNode 4-nNode 5-nNode 1 2 3 4 5 6 7 8 9];
        nFile = 2;
     case 421
        nNode = 1;%�ڵ���
        loadID = [1-nNode 1 2];
        nFile = 1;
     case 422
        nNode = 1;%�ڵ���
        loadID = [1-nNode 1];
        nFile = 1;
     case 423
        nNode = 5;%�ڵ���
        loadID = [1-nNode 2-nNode 3-nNode 4-nNode 5-nNode 1 2 3 4 5 6 7 8 9];
        nFile = 2;
     case 424
        nNode = 5;%�ڵ���
        loadID = [1-nNode 2-nNode 3-nNode 4-nNode 5-nNode 1 2 3 4 5 6 7 8 9];
        nFile = 2;
     case 425
        nNode = 5;%�ڵ��������ز���ģ��
        loadID = [1-nNode 2-nNode 3-nNode 4-nNode 5-nNode 1 2 3 4 5 6 7 8 9];
        nFile = 2;
     case 426
        nNode = 1;%�ڵ��������ز�������
        loadID = [1-nNode];
        nFile = 1;
     case 427
        nNode = 1;%�ڵ��������ز�������ϵͳ��������
        loadID = [1-nNode];
        nFile = 1;
     case 428  
        nNode = 33;%�ڵ�������΢������ģ��
        loadID = [1-nNode 2-nNode 3-nNode];
        nFile = 10;
     case 429   
        nNode = 34;%�ڵ�������΢������ģ��
        loadID = [1-nNode 2-nNode 3-nNode];
        nFile = 10;
     case 430   
        nNode = 1;%�ڵ�������������ģ���������
        loadID = [1-nNode];
        nFile = 1;
     case 431
        nNode = 1;%�ڵ�������΢����������ϵͳ��������
        loadID = [1-nNode];
        nFile = 1;
     case 432  
        nNode = 62;%�ڵ�������΢������ģ��
        loadID = [1-nNode 2-nNode 3-nNode];
        nFile = 16;
     case 433  
        nNode = 54;%�ڵ�����ȥ����ѹ���ļ�΢������ģ��
        loadID = [1-nNode 2-nNode 3-nNode];
        nFile = 14;   
     case 434  
        nNode = 54;%�ڵ�����ȥ����ѹ���ļ�΢������ģ�ͣ����߶�·���г����ߺ͹��1�ı�������
        loadID = [1-nNode 2-nNode 3-nNode];
        nFile = 14; 
    case 435  
        nNode = 54;%�ڵ�����ȥ����ѹ���ļ�΢������ģ�ͣ�PCC��·�������ı�������
        loadID = [1-nNode 2-nNode 3-nNode];
        nFile = 14; 
end
%% ��ȡPSCAD����
result_PSC = [];
for k = 1:nFile
    if k<9.5
        stg = sprintf('./result/pscadData/case%d_0%d.out',iCase,k);
    else
        stg = sprintf('./result/pscadData/case%d_%d.out',iCase,k);
    end
    result_temp = load(stg);
    t_PSC=result_temp(:,1);
    dT = t_PSC(2)-t_PSC(1);
    result_PSC = [result_PSC result_temp(:,2:end)*1000];
end

%��ȡCPP����
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t_CPP=result_temp(:,1);%ʱ������
result_CPP=result_temp(:,2:end);%CPP����

%%

%��ͼ����
if plotFlag==1
   index = [loadID+nNode];%��Ҫ��ͼ���б��,�ֱ�Ϊ���ص����������������Դ
    ss = size(index);
    for k = 1:ss(2)
        figure(index(k));
        subplot(2,1,1);
        plot(t_CPP,result_CPP(:,index(k)),'--',t_PSC,result_PSC(:,index(k)),':');
        ymax=max(max(result_PSC(:,index(k))),max(result_CPP(:,index(k))));
        ymin=min(min(result_PSC(:,index(k))),min(result_CPP(:,index(k))));
        yh=ymax-ymin;
        ymax=ymax+yh/2;
        ymin=ymin-yh/2;
        ylim([ymin ymax]);
        legend('CPP result','PSC result')
        subplot(2,1,2);
        plot(t_CPP,result_CPP(:,index(k))-result_PSC((t_CPP(1)-t_PSC(1))/dT+1:end,index(k)));
    end
end

%%
% ����ָ��
% [VsPSCMag VsPSCPha]=cal_mag_theta(result_PSC(end,1:3)');
% [IsPSCMag IsPSCPha]=cal_mag_theta(result_PSC(end,7:9)');
% [VrPSCMag VrPSCPha]=cal_mag_theta(result_PSC(end,4:6)');
% [VsCPPMag VsCPPPha]=cal_mag_theta(result_CPP(end,1:3)');
% [IsCPPMag IsCPPPha]=cal_mag_theta(result_CPP(end,7:9)');
% [VrCPPMag VrCPPPha]=cal_mag_theta(result_CPP(end,4:6)');
% 
% VsErrMag=abs(VsPSCMag-VsCPPMag)/VsPSCMag*100;
% VrErrMag=abs(VrPSCMag-VrCPPMag)/VrPSCMag*100;
% IsErrMag=abs(IsPSCMag-IsCPPMag)/IsPSCMag*100;
% VsErrPha=(VsPSCPha-VsCPPPha)*180/pi;
% VrErrPha=(VrPSCPha-VrCPPPha)*180/pi;
% IsErrPha=(IsPSCPha-IsCPPPha)*180/pi;









