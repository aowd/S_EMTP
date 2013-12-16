clc
clear all;
close all;
format long
stg=strcat('================',datestr(now),'================');
disp(stg);
%% 算例信息
% iCase = 76;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 8;%节点数
% loadID = [1 13 4-nNode 7-nNode];
% nFile = 3;

% iCase = 7636;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 32;%节点数
% loadID = [56 1 7-nNode 32-nNode];
% nFile = 10;

% iCase = 1103;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 12;%节点数
% loadID = [1 22 1-nNode 2-nNode];
% nFile = 5;

% iCase = 1110;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 44;%节点数
% loadID = [107 120 3-nNode 4-nNode];
% nFile = 17;

% iCase = 71103;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 20;%节点数
% loadID = [45 1 50 7-nNode 8-nNode];
% nFile = 8;

% iCase = 71104;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 25;%节点数
% loadID = [45 1 55 7-nNode 8-nNode];
% nFile = 9;

% iCase = 71110;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 68;%节点数
% loadID = [160 45 168 7-nNode 32-nNode 35-nNode];
% nFile = 25;

% iCase = 71111;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 83;%节点数
% loadID = [174 2 183 7-nNode 32-nNode];
% nFile = 28;

% iCase = 7604;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 8;%节点数
% loadID = 13;%负载电阻支路编号
% nFile = 3;

% iCase = 11111;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = 4;%负载电阻支路编号
% nFile = 1;

% iCase = 12003;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 24;%节点数
% loadID = [1 1-nNode];
% nFile = 6;

% iCase = 120071;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 13;%节点数
% loadID = [2 13-nNode];
% nFile = 4;

% iCase = 12007;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 25;%节点数
% loadID = [1 2 25-nNode];
% nFile = 7;

% iCase = 12023;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 49;%节点数
% loadID = [1 80 1-nNode];
% nFile = 13;

% iCase = 15001;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 30;%节点数
% loadID = [1 1-nNode];
% nFile = 9;

% iCase = 150011;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 60;%节点数
% loadID = [1 1-nNode];
% nFile = 17;

% iCase = 15002;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 45;%节点数
% loadID = [1 1-nNode];
% nFile = 12;

% iCase = 15003;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 31;%节点数
% loadID = [1 1-nNode 31-nNode];
% nFile = 10;

% iCase = 15004;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 46;%节点数
% loadID = [1 1-nNode 46-nNode];
% nFile = 13;

% iCase = 15005;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 45;%节点数
% loadID = [1 1-nNode 31-nNode];
% nFile = 12;

% iCase = 15006;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 60;%节点数
% loadID = [1 1-nNode];
% nFile = 15;

% iCase = 15007;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 38;%节点数
% loadID = [1 109 1-nNode 31-nNode];
% nFile = 16;

% iCase = 150071;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 68;%节点数
% loadID = [1 1-nNode];
% nFile = 22;

% iCase = 666666;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 56;%节点数
% loadID = [1 109 156];
% nFile = 22;

% iCase = 603;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 6;%节点数
% loadID = [1 1-nNode];
% nFile = 2;

% iCase = 101;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 60;%节点数
% loadID = [1-nNode 1];
% nFile = 17;

% iCase = 102;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 24;%节点数
% loadID = [1-nNode 1];
% nFile = 6;

% iCase = 103;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 33;%节点数
% loadID = [1-nNode 1];
% nFile = 10;

% iCase = 104;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 68;%节点数
% loadID = [1-nNode 1];
% nFile = 22;

% iCase = 105;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 68;%节点数
% loadID = [1-nNode 1];
% nFile = 22;

% iCase = 106;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 49;%节点数
% loadID = [1-nNode 1];
% nFile = 13;

% iCase = 107;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 83;%节点数
% loadID = [1-nNode 1];
% nFile = 28;

% iCase = 108;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 98;%节点数
% loadID = [1-nNode 67-nNode 87-nNode 1 109 198];
% nFile = 31;

% iCase = 109;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 80;%节点数
% loadID = [1-nNode 37-nNode 69-nNode 1 109 180];
% nFile = 27;

% iCase = 132;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 1;%节点数
% loadID = [1-nNode 1 2];
% nFile = 1;

% iCase = 201;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 6;%节点数
% loadID = [1-nNode 1];
% nFile = 2;

% iCase = 202;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 9;%节点数
% loadID = [1-nNode 7-nNode 1 15];
% nFile = 3;

% iCase = 204;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [1-nNode 1];
% nFile = 1;

% iCase = 205;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [1-nNode 1];
% nFile = 1;

% iCase = 206;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 9;%节点数
% loadID = [4-nNode 5];
% nFile = 3;

% iCase = 207;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [1-nNode 1];
% nFile = 1;

% iCase = 208;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 9;%节点数
% loadID = [4-nNode 5];
% nFile = 3;

% iCase = 209;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 6;%节点数
% loadID = [4-nNode 5];
% nFile = 2;

% iCase = 210;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 15;%节点数
% loadID = [5-nNode 2 15];
% nFile = 4;

% iCase = 211;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 15;%节点数
% loadID = [5-nNode 1 6];
% nFile = 4;

% iCase = 212;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 6;%节点数
% loadID = [4-nNode 2 5];
% nFile = 2;

% iCase = 213;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 15;%节点数
% loadID = [5-nNode 15];
% nFile = 4;

% iCase = 214;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 18;%节点数
% loadID = [5-nNode 15 1 24];
% nFile = 5;

% iCase = 216;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 15;%节点数
% loadID = [2-nNode 2];
% nFile = 4;
% 
% iCase = 219;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [1-nNode 2-nNode 3-nNode 3];
% nFile = 1;
% 
% iCase = 220;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 2;%节点数
% loadID = [2-nNode];
% nFile = 1;

% iCase = 220;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 2;%节点数
% loadID = [2-nNode 2];
% nFile = 1;

% iCase = 221;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 2;%节点数
% loadID = [1-nNode 2-nNode 1 2 3];
% nFile = 1;

% iCase = 222;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 2;%节点数
% loadID = [1-nNode 2-nNode 1 2 3];
% nFile = 1;

% iCase = 225;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [2-nNode 1 2 3];
% nFile = 1;

% iCase = 227;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [1-nNode 2-nNode 2 3 4];
% nFile = 1;


% iCase = 228;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 4;%节点数
% loadID = [1-nNode 2-nNode 2 3 4];
% nFile = 2;

% iCase = 230;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 4;%节点数
% loadID = [2-nNode 3-nNode 1 2 3 10];
% nFile = 2;

% iCase = 231;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 4;%节点数
% loadID = [1-nNode 2-nNode 3-nNode 4-nNode 2 3 5 8 10];
% nFile = 2;

% iCase = 232;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 11;%节点数
% loadID = [1-nNode 2-nNode 3-nNode 4-nNode 1 15];
% nFile = 4;

% iCase = 243;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 1;%节点数
% loadID = [1-nNode];
% nFile = 1;

% iCase = 243;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 1;%节点数
% loadID = [1-nNode];
% nFile = 1;
% 
% iCase = 245;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 6;%节点数
% loadID = [1-nNode 3 4 5 6 7 8 9 10 13];
% nFile = 2;

% iCase = 248;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 3;%节点数
% loadID = [3-nNode 2];
% nFile = 1;

% iCase = 257;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 1;%节点数
% loadID = [1-nNode];
% nFile = 1;

% iCase = 259;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 1;%节点数
% loadID = [1-nNode];
% nFile = 1;

% iCase = 260;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 18;%节点数
% loadID = [1-nNode 1 2 4 6];
% nFile = 5;


% iCase = 261;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 23;%节点数
% loadID = [1-nNode 4-nNode 14 17];
% nFile = 8;

% iCase = 262;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 11;%节点数
% loadID = [1-nNode 4-nNode 14 17];
% nFile = 4;

% iCase = 263;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 1;%节点数
% loadID = [1-nNode];
% nFile = 1;

% iCase = 264;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 14;%节点数
% loadID = [1-nNode 4-nNode 1 7 14 17];
% nFile = 5;

% iCase = 270;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 11;%节点数
% loadID = [1-nNode 4 5 16];
% nFile = 4;

iCase = 271;%算例编号
plotFlag = 1;%作图标志位
nNode = 5;%节点数
loadID = [1 5];
nFile = 2;
deltaT_ratio = 1;

% iCase = 258;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 11;%节点数
% loadID = [1 19];
% nFile = 4;
% deltaT_ratio = 4;

% iCase = 272;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 5;%节点数
% loadID = [1 5];
% nFile = 2;
% deltaT_ratio = 1;

% iCase = 273;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 8;%节点数
% loadID = [4-nNode 6-nNode 1];
% nFile = 3;
% deltaT_ratio = 1;

% iCase = 500;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 10;%节点数
% loadID = [1-nNode 1];
% nFile = 3;

% iCase = 502;%算例编号
% plotFlag = 1;%作图标志位
% nNode = 12;%节点数
% loadID = [1-nNode 1];
% nFile = 5;


%% 读取PSCAD数据
result_PSC = [];
for k = 1:nFile
    if k<9.5
        stg = sprintf('./result/pscadData/case%d_0%d.out',iCase,k);
    else
        stg = sprintf('./result/pscadData/case%d_%d.out',iCase,k);
    end
    result_temp = load(stg);
    time_PSC = result_temp(:,1);
    result_PSC = [result_PSC result_temp(:,2:end)*1000];
end

%读取CPP数据
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%时间序列
result_CPP=result_temp(:,2:end);%CPP数据

% result_PSC = result_PSC(1:size(result_CPP,1),:);
idx = 1;
t_tmp = time_PSC(idx);
while t_tmp ~= t(1)
    idx = idx + 1;
    t_tmp = time_PSC(idx);
end
result_PSC = result_PSC(idx:deltaT_ratio:idx+size(result_CPP,1)*deltaT_ratio-1,:);

% %数据对比
% temp=abs(result_CPP-result_PSC);
% err=max(max(temp));
% [i,j]=find(temp==err);
% stg=sprintf('绝对误差的最大值为：err = %g',err);
% disp(stg);
% if j<nNode+0.5
%     stg=sprintf('最大误差来自节点%d的电压，出现于counter=%d，time=%ds的时刻。',j(1),i(1),(i(1)-1)*50e-6);
% else
%     stg=sprintf('最大误差来自支路%d的电流，出现于counter=%d，time=%ds的时刻。',j(1)-nNode,i(1),(i(1)-1)*50e-6);
% end
% disp(stg);
% 
% stg=strcat('====================================================');
% disp(stg);

%作图环节
if plotFlag==1
%     index = [loadID+nNode j(1)];%需要作图的列编号,分别为负载电阻电流和最大误差来源
   index = [loadID+nNode];%需要作图的列编号,分别为负载电阻电流和最大误差来源
    ss = size(index);
    for k = 1:ss(2)
        figure(index(k));
        subplot(2,1,1)
        plot(t,result_PSC(:,index(k)),':',t,result_CPP(:,index(k)),'--');
%         ymax=max(max(result_PSC(:,index(k))),max(result_CPP(:,index(k))));
%         ymin=min(min(result_PSC(:,index(k))),min(result_CPP(:,index(k))));
%         yh=ymax-ymin;
%         ymax=ymax+yh/2;
%         ymin=ymin-yh/2;
%         ylim([ymin ymax]);
        xlim([t(1) t(end)]);
        legend('PSC result','CPP result')
        subplot(2,1,2)
        plot(t,result_CPP(:,index(k))-result_PSC(:,index(k)));
        xlim([t(1) t(end)]);
    end
end


%% debug
% CPP = result_CPP(:,1)-result_CPP(:,2);
% PSC = result_PSC(:,1)-result_PSC(:,2);
% 
% figure(1000);
% subplot(2,1,1)
% plot(t,CPP,'--',t,PSC,':');
% ymax=max(max(PSC),max(CPP));
% ymin=min(min(PSC),min(CPP));
% yh=ymax-ymin;
% ymax=ymax+yh/2;
% ymin=ymin-yh/2;
% ylim([ymin ymax]);
% legend('CPP result','PSC result')
% subplot(2,1,2)
% plot(t,CPP-PSC);
% 
% 
% CPP = result_CPP(:,49)-result_CPP(:,56);
% PSC = result_PSC(:,49)-result_PSC(:,56);
% 
% figure(2000);
% subplot(2,1,1)
% plot(t,CPP,'--',t,PSC,':');
% ymax=max(max(PSC),max(CPP));
% ymin=min(min(PSC),min(CPP));
% yh=ymax-ymin;
% ymax=ymax+yh/2;
% ymin=ymin-yh/2;
% ylim([ymin ymax]);
% legend('CPP result','PSC result')
% subplot(2,1,2)
% plot(t,CPP-PSC);










