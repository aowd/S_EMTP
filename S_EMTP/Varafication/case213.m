clear all;
clc;

%% ��ȡ�ļ�
% ������Ϣ
iCase = 213;
nFile = 4;
deltaT = 50;
deltaT_Init = 50;

% ��ȡPSCAD�����ļ�
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

% ��ȡCPP�����ļ�
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%ʱ������
result_CPP=result_temp(:,2:end);%CPP����

% ����PSCAD���ݣ�ʹ֮��CPP���ݶ�Ӧ
steps = deltaT/deltaT_Init;
idx = 1;
t_tmp = time_PSC(idx);
while t_tmp ~= t(1)
    idx = idx + 1;
    t_tmp = time_PSC(idx);
end
result_PSC = result_PSC(idx:steps:idx+size(result_CPP,1)*steps-steps,:);

% ����ָ��������ʼ��
testStart = int32(size(t,1)*0.6);

%% ѡ����Ҫ�۲�Ĳ���
% ���ӵ���
% Iabc_PSC = result_PSC(:,16:18);
% Iabc_CPP = result_CPP(:,16:18);

% ת�ӵ���
Iabc_PSC = result_PSC(:,19:21);
Iabc_CPP = result_CPP(:,19:21);

%% �˲�
% �˲�������ϵ��
load Num_ac_50us.mat
Num_ac = Num_ac_50us;

% �����е���abc˲ʱֵ
Iabc_PSC_filter = filter(Num_ac,1,Iabc_PSC);
Iabc_CPP_filter = filter(Num_ac,1,Iabc_CPP);

% ��ֵ��Ǽ���
for k = 1:size(Iabc_PSC_filter,1)
    [Mag_PSC_filter(k) phase_PSC_filter(k)] = cal_mag_theta(Iabc_PSC_filter(k,:)');
    [Mag_CPP_filter(k) phase_CPP_filter(k)] = cal_mag_theta(Iabc_CPP_filter(k,:)');
end

delta_Mag = abs( ( Mag_CPP_filter(testStart:end) - Mag_PSC_filter(testStart:end) ) ./ Mag_PSC_filter(testStart:end) );
Err_Mag = max(delta_Mag)*100

% ��ǲ����
for k = 1:size(Iabc_PSC_filter,1)
    delta_phase(k) = phase_CPP_filter(k) - phase_PSC_filter(k);
    if delta_phase(k) > pi
        delta_phase(k) = delta_phase(k) - 2*pi;
    end
    if delta_phase(k) < -pi
        delta_phase(k) = delta_phase(k) + 2*pi;
    end
    delta_phase(k) = delta_phase(k)/pi*180;
end

Err_phase = max(abs(delta_phase(testStart:end)))

% ʵ���鲿����
Ir_PSC_filter = (2/3)*Iabc_PSC_filter(:,1) - (1/3)*Iabc_PSC_filter(:,2) ...
     - (1/3)*Iabc_PSC_filter(:,3);
Ii_PSC_filter = (1/sqrt(3))*(Iabc_PSC_filter(:,2)-Iabc_PSC_filter(:,3));

Ir_CPP_filter = (2/3)*Iabc_CPP_filter(:,1) - (1/3)*Iabc_CPP_filter(:,2) ...
     - (1/3)*Iabc_CPP_filter(:,3);
Ii_CPP_filter = (1/sqrt(3))*(Iabc_CPP_filter(:,2)-Iabc_CPP_filter(:,3));

TVE = sqrt( ( (Ir_PSC_filter-Ir_CPP_filter).^2+(Ii_PSC_filter-Ii_CPP_filter).^2 ) ./ ( Ir_PSC_filter.^2 + Ii_PSC_filter.^2 ) );
Err_iabc = max(TVE(testStart:end))*100

% 0���������
i0_PSC = (1/3)*(Iabc_PSC_filter(:,1)+Iabc_PSC_filter(:,2)+Iabc_PSC_filter(:,3));
i0_CPP = (1/3)*(Iabc_CPP_filter(:,1)+Iabc_CPP_filter(:,2)+Iabc_CPP_filter(:,3));
delta_i0 = abs( ( i0_CPP(testStart:end) - i0_PSC(testStart:end) ) ./ i0_PSC(testStart:end) );
Err_i0 = max(delta_i0)*100

%% ��ͼ�Ƚ�
figure(1)
plot(t,Iabc_PSC_filter,t,Iabc_CPP_filter);
 
figure(2)
plot(t,Mag_PSC_filter,t,Mag_CPP_filter);

figure(3)
plot(t,delta_phase);

figure(4)
plot(t,i0_PSC,t,i0_CPP)
