clc
clear all;
close all;
format long
stg=strcat('================',datestr(now),'================');
disp(stg);
%% ������Ϣ
% iCase = 76;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 8;%�ڵ���
% loadID = [1 13 4-nNode 7-nNode];
% nFile = 3;

% iCase = 7636;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 32;%�ڵ���
% loadID = [56 1 7-nNode 32-nNode];
% nFile = 10;

% iCase = 1103;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 12;%�ڵ���
% loadID = [1 22 1-nNode 2-nNode];
% nFile = 5;

% iCase = 1110;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 44;%�ڵ���
% loadID = [107 120 3-nNode 4-nNode];
% nFile = 17;

% iCase = 71103;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 20;%�ڵ���
% loadID = [45 1 50 7-nNode 8-nNode];
% nFile = 8;

% iCase = 71104;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 25;%�ڵ���
% loadID = [45 1 55 7-nNode 8-nNode];
% nFile = 9;

% iCase = 71110;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 68;%�ڵ���
% loadID = [160 45 168 7-nNode 32-nNode 35-nNode];
% nFile = 25;

% iCase = 71111;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 83;%�ڵ���
% loadID = [174 2 183 7-nNode 32-nNode];
% nFile = 28;

% iCase = 7604;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 8;%�ڵ���
% loadID = 13;%���ص���֧·���
% nFile = 3;

% iCase = 11111;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = 4;%���ص���֧·���
% nFile = 1;

% iCase = 12003;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 24;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 6;

% iCase = 120071;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 13;%�ڵ���
% loadID = [2 13-nNode];
% nFile = 4;

% iCase = 12007;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 25;%�ڵ���
% loadID = [1 2 25-nNode];
% nFile = 7;

% iCase = 12023;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 49;%�ڵ���
% loadID = [1 80 1-nNode];
% nFile = 13;

% iCase = 15001;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 30;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 9;

% iCase = 150011;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 60;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 17;

% iCase = 15002;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 45;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 12;

% iCase = 15003;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 31;%�ڵ���
% loadID = [1 1-nNode 31-nNode];
% nFile = 10;

% iCase = 15004;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 46;%�ڵ���
% loadID = [1 1-nNode 46-nNode];
% nFile = 13;

% iCase = 15005;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 45;%�ڵ���
% loadID = [1 1-nNode 31-nNode];
% nFile = 12;

% iCase = 15006;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 60;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 15;

% iCase = 15007;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 38;%�ڵ���
% loadID = [1 109 1-nNode 31-nNode];
% nFile = 16;

% iCase = 150071;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 68;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 22;

% iCase = 666666;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 56;%�ڵ���
% loadID = [1 109 156];
% nFile = 22;

% iCase = 603;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 6;%�ڵ���
% loadID = [1 1-nNode];
% nFile = 2;

% iCase = 101;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 60;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 17;

% iCase = 102;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 24;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 6;

% iCase = 103;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 33;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 10;

% iCase = 104;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 68;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 22;

% iCase = 105;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 68;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 22;

% iCase = 106;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 49;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 13;

% iCase = 107;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 83;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 28;

% iCase = 108;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 98;%�ڵ���
% loadID = [1-nNode 67-nNode 87-nNode 1 109 198];
% nFile = 31;

% iCase = 109;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 80;%�ڵ���
% loadID = [1-nNode 37-nNode 69-nNode 1 109 180];
% nFile = 27;

% iCase = 132;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 1;%�ڵ���
% loadID = [1-nNode 1 2];
% nFile = 1;

% iCase = 201;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 6;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 2;

% iCase = 202;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 9;%�ڵ���
% loadID = [1-nNode 7-nNode 1 15];
% nFile = 3;

% iCase = 204;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 1;

% iCase = 205;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 1;

% iCase = 206;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 9;%�ڵ���
% loadID = [4-nNode 5];
% nFile = 3;

% iCase = 207;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 1;

% iCase = 208;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 9;%�ڵ���
% loadID = [4-nNode 5];
% nFile = 3;

% iCase = 209;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 6;%�ڵ���
% loadID = [4-nNode 5];
% nFile = 2;

% iCase = 210;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 15;%�ڵ���
% loadID = [5-nNode 2 15];
% nFile = 4;

% iCase = 211;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 15;%�ڵ���
% loadID = [5-nNode 1 6];
% nFile = 4;

% iCase = 212;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 6;%�ڵ���
% loadID = [4-nNode 2 5];
% nFile = 2;

% iCase = 213;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 15;%�ڵ���
% loadID = [5-nNode 15];
% nFile = 4;

% iCase = 214;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 18;%�ڵ���
% loadID = [5-nNode 15 1 24];
% nFile = 5;

% iCase = 216;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 15;%�ڵ���
% loadID = [2-nNode 2];
% nFile = 4;
% 
% iCase = 219;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [1-nNode 2-nNode 3-nNode 3];
% nFile = 1;
% 
% iCase = 220;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 2;%�ڵ���
% loadID = [2-nNode];
% nFile = 1;

% iCase = 220;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 2;%�ڵ���
% loadID = [2-nNode 2];
% nFile = 1;

% iCase = 221;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 2;%�ڵ���
% loadID = [1-nNode 2-nNode 1 2 3];
% nFile = 1;

% iCase = 222;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 2;%�ڵ���
% loadID = [1-nNode 2-nNode 1 2 3];
% nFile = 1;

% iCase = 225;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [2-nNode 1 2 3];
% nFile = 1;

% iCase = 227;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [1-nNode 2-nNode 2 3 4];
% nFile = 1;


% iCase = 228;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 4;%�ڵ���
% loadID = [1-nNode 2-nNode 2 3 4];
% nFile = 2;

% iCase = 230;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 4;%�ڵ���
% loadID = [2-nNode 3-nNode 1 2 3 10];
% nFile = 2;

% iCase = 231;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 4;%�ڵ���
% loadID = [1-nNode 2-nNode 3-nNode 4-nNode 2 3 5 8 10];
% nFile = 2;

% iCase = 232;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 11;%�ڵ���
% loadID = [1-nNode 2-nNode 3-nNode 4-nNode 1 15];
% nFile = 4;

% iCase = 243;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 1;%�ڵ���
% loadID = [1-nNode];
% nFile = 1;

% iCase = 243;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 1;%�ڵ���
% loadID = [1-nNode];
% nFile = 1;
% 
% iCase = 245;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 6;%�ڵ���
% loadID = [1-nNode 3 4 5 6 7 8 9 10 13];
% nFile = 2;

% iCase = 248;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 3;%�ڵ���
% loadID = [3-nNode 2];
% nFile = 1;

% iCase = 257;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 1;%�ڵ���
% loadID = [1-nNode];
% nFile = 1;

% iCase = 259;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 1;%�ڵ���
% loadID = [1-nNode];
% nFile = 1;

% iCase = 260;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 18;%�ڵ���
% loadID = [1-nNode 1 2 4 6];
% nFile = 5;


% iCase = 261;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 23;%�ڵ���
% loadID = [1-nNode 4-nNode 14 17];
% nFile = 8;

% iCase = 262;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 11;%�ڵ���
% loadID = [1-nNode 4-nNode 14 17];
% nFile = 4;

% iCase = 263;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 1;%�ڵ���
% loadID = [1-nNode];
% nFile = 1;

% iCase = 264;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 14;%�ڵ���
% loadID = [1-nNode 4-nNode 1 7 14 17];
% nFile = 5;

% iCase = 270;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 11;%�ڵ���
% loadID = [1-nNode 4 5 16];
% nFile = 4;

iCase = 271;%�������
plotFlag = 1;%��ͼ��־λ
nNode = 5;%�ڵ���
loadID = [1 5];
nFile = 2;
deltaT_ratio = 1;

% iCase = 258;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 11;%�ڵ���
% loadID = [1 19];
% nFile = 4;
% deltaT_ratio = 4;

% iCase = 272;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 5;%�ڵ���
% loadID = [1 5];
% nFile = 2;
% deltaT_ratio = 1;

% iCase = 273;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 8;%�ڵ���
% loadID = [4-nNode 6-nNode 1];
% nFile = 3;
% deltaT_ratio = 1;

% iCase = 500;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 10;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 3;

% iCase = 502;%�������
% plotFlag = 1;%��ͼ��־λ
% nNode = 12;%�ڵ���
% loadID = [1-nNode 1];
% nFile = 5;


%% ��ȡPSCAD����
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

%��ȡCPP����
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%ʱ������
result_CPP=result_temp(:,2:end);%CPP����

% result_PSC = result_PSC(1:size(result_CPP,1),:);
idx = 1;
t_tmp = time_PSC(idx);
while t_tmp ~= t(1)
    idx = idx + 1;
    t_tmp = time_PSC(idx);
end
result_PSC = result_PSC(idx:deltaT_ratio:idx+size(result_CPP,1)*deltaT_ratio-1,:);

% %���ݶԱ�
% temp=abs(result_CPP-result_PSC);
% err=max(max(temp));
% [i,j]=find(temp==err);
% stg=sprintf('�����������ֵΪ��err = %g',err);
% disp(stg);
% if j<nNode+0.5
%     stg=sprintf('���������Խڵ�%d�ĵ�ѹ��������counter=%d��time=%ds��ʱ�̡�',j(1),i(1),(i(1)-1)*50e-6);
% else
%     stg=sprintf('����������֧·%d�ĵ�����������counter=%d��time=%ds��ʱ�̡�',j(1)-nNode,i(1),(i(1)-1)*50e-6);
% end
% disp(stg);
% 
% stg=strcat('====================================================');
% disp(stg);

%��ͼ����
if plotFlag==1
%     index = [loadID+nNode j(1)];%��Ҫ��ͼ���б��,�ֱ�Ϊ���ص����������������Դ
   index = [loadID+nNode];%��Ҫ��ͼ���б��,�ֱ�Ϊ���ص����������������Դ
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










