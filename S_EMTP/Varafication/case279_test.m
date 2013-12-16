clear all;
clc;

%% 读取文件
% 算例信息
iCase = 279;
nFile = 3;
deltaT = 5;

% 读取CPP数据文件
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%时间序列
result_CPP=result_temp(:,2:end);%CPP数据

%% 选择需要观察的波形
Iabc_CPP = result_CPP(:,12:14);

%% 滤波
% 滤波器分子系数
load Num_ac_5us.mat
Num_ac = Num_ac_5us;

% 三相电感电流abc瞬时值
Iabc_CPP_filter = filter(Num_ac,1,Iabc_CPP);

i0_CPP = (1/3)*(Iabc_CPP_filter(:,1)+Iabc_CPP_filter(:,2)+Iabc_CPP_filter(:,3));

%% 作图比较
figure(1)
plot(t,i0_CPP)

