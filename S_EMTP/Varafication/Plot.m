clc
clear all;
close all;
format short

% iCase = 303;%算例编号
% nNode = 130;%节点数
% VoltageNum = [1 49 100];
% CurrentNum = [300];

% iCase = 205;%算例编号
% nNode = 3;%节点数
% VoltageNum = [1];
% CurrentNum = [1];
% 
% iCase = 207;%算例编号
% nNode = 3;%节点数
% VoltageNum = [1];
% CurrentNum = [1];

% iCase = 216;%算例编号
% nNode = 15;%节点数
% VoltageNum = [5];
% CurrentNum = [14];

iCase = 221;%算例编号
nNode = 2;%节点数
VoltageNum = [2];
CurrentNum = [1 2 3 4];


stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%时间序列
result_CPP=result_temp(:,2:end);%CPP数据

index = [VoltageNum CurrentNum+nNode];%需要作图的列编号,分别为负载电阻电流和最大误差来源
ss = size(index);
for k = 1:ss(2)
    figure(index(k));
    plot(t,result_CPP(:,index(k)));
    ymax=max(result_CPP(:,index(k)));
    ymin=min(result_CPP(:,index(k)));
    yh=ymax-ymin;
    ymax=ymax+yh/2;
    ymin=ymin-yh/2;
    ylim([ymin ymax]);
    legend('CPP result')
end