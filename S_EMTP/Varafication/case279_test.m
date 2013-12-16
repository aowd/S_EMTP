clear all;
clc;

%% ��ȡ�ļ�
% ������Ϣ
iCase = 279;
nFile = 3;
deltaT = 5;

% ��ȡCPP�����ļ�
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%ʱ������
result_CPP=result_temp(:,2:end);%CPP����

%% ѡ����Ҫ�۲�Ĳ���
Iabc_CPP = result_CPP(:,12:14);

%% �˲�
% �˲�������ϵ��
load Num_ac_5us.mat
Num_ac = Num_ac_5us;

% �����е���abc˲ʱֵ
Iabc_CPP_filter = filter(Num_ac,1,Iabc_CPP);

i0_CPP = (1/3)*(Iabc_CPP_filter(:,1)+Iabc_CPP_filter(:,2)+Iabc_CPP_filter(:,3));

%% ��ͼ�Ƚ�
figure(1)
plot(t,i0_CPP)

