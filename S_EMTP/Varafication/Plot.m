clc
clear all;
close all;
format short

% iCase = 303;%�������
% nNode = 130;%�ڵ���
% VoltageNum = [1 49 100];
% CurrentNum = [300];

% iCase = 205;%�������
% nNode = 3;%�ڵ���
% VoltageNum = [1];
% CurrentNum = [1];
% 
% iCase = 207;%�������
% nNode = 3;%�ڵ���
% VoltageNum = [1];
% CurrentNum = [1];

% iCase = 216;%�������
% nNode = 15;%�ڵ���
% VoltageNum = [5];
% CurrentNum = [14];

iCase = 221;%�������
nNode = 2;%�ڵ���
VoltageNum = [2];
CurrentNum = [1 2 3 4];


stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%ʱ������
result_CPP=result_temp(:,2:end);%CPP����

index = [VoltageNum CurrentNum+nNode];%��Ҫ��ͼ���б��,�ֱ�Ϊ���ص����������������Դ
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