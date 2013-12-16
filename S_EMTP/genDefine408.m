% =============================================================
%         ��case 408Ϊ����������һϵ�е����ķ������
%                 ����֮�䲢������
% =============================================================
clear all
clc

% �����ʼ�ļ�
load ./caseDfn/case408.txt
load ./caseDfn/case408ctrl.txt
load ./caseDfn/case408msr.txt

% ����ڵ���
nElecNode = 17;
nCtrlNode = 193;
nElecBranch = size(case408,1);

% ����DFIG����
nDFIG = 72;

%% ���ɵ���ϵͳ�������
elecMatrix = case408;
temp = case408;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,2) ~= 0
            temp(i,2) = temp(i,2) + nElecNode;
        end
        if temp(i,3) ~= 0
            temp(i,3) = temp(i,3) + nElecNode;
        end
        if temp(i,1)==24
            temp(i,7)=temp(i,7) + nCtrlNode;
        end
    end
    elecMatrix = [elecMatrix; temp];
end

%% ��������ϵͳ�������
msrMatrix = case408msr;
temp = case408msr;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,1) == 1
            temp(i,2) = temp(i,2) + nElecNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
        if temp(i,1) == 3 || temp(i,1) == 4 ||temp(i,1) == 5 
            temp(i,2) = temp(i,2) + nElecBranch;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
    end
    msrMatrix = [msrMatrix; temp];
end

%% ���ɿ���ϵͳ�������
ctrlMatrix = case408ctrl;
temp = case408ctrl;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,1) == 3 || temp(i,1) == 5
            temp(i,2) = temp(i,2) + nCtrlNode;
        elseif temp(i,1) == 2 || temp(i,1) == 7 || temp(i,1) == 14
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        elseif temp(i,1) == 1
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
        elseif temp(i,1) == 12 || temp(i,1) == 13
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
            temp(i,5) = temp(i,5) + nCtrlNode;
        elseif temp(i,1) == 8 || temp(i,1) == 9 || temp(i,1) == 10 || temp(i,1) == 11
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
            temp(i,5) = temp(i,5) + nCtrlNode;
            temp(i,6) = temp(i,6) + nCtrlNode;
        end
    end
    ctrlMatrix = [ctrlMatrix; temp];
end

%% �洢���
stg = strcat(['./caseDfn/case408',num2str(nDFIG)],'.txt');
save(stg,'-ascii','elecMatrix');
stg = strcat(['./caseDfn/case408',num2str(nDFIG)],'ctrl.txt');
save(stg,'-ascii','ctrlMatrix');
stg = strcat(['./caseDfn/case408',num2str(nDFIG)],'msr.txt');
save(stg,'-ascii','msrMatrix');