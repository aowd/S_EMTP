% clear all;
% clc;

%% 读取文件
% 算例信息
iCase = 273;
nFile = 3;
deltaT = 5;
testStart = 10000;

% 读取PSCAD数据文件
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

% 读取CPP数据文件
stg = sprintf('./result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%时间序列
result_CPP=result_temp(:,2:end);%CPP数据

% 处理PSCAD数据，使之与CPP数据对应
steps = deltaT/5;
idx = 1;
t_tmp = time_PSC(idx);
while t_tmp ~= t(1)
    idx = idx + 1;
    t_tmp = time_PSC(idx);
end
result_PSC = result_PSC(idx:steps:idx+size(result_CPP,1)*steps-steps,:);

%% 选择需要观察的波形
Udc_PSC = result_PSC(:,4)-result_PSC(:,5);
Udc_CPP = result_CPP(:,4)-result_CPP(:,5);

Iabc_PSC = result_PSC(:,12:14);
Iabc_CPP = result_CPP(:,12:14);

Iabc_PSC = result_PSC(:,17:19);
Iabc_CPP = result_CPP(:,17:19);

%% 滤波
% 滤波器分子系数
% load Num_ac.mat

% 直流电压
Udc_PSC_filter = filter(Num_ac,1,Udc_PSC);
Udc_CPP_filter = filter(Num_ac,1,Udc_CPP);

% delta_Udc = abs( ( Udc_CPP_filter(10000:end) - Udc_PSC_filter(10000:end) ) ./ Udc_PSC_filter(10000:end) );
delta_Udc = abs( ( Udc_CPP_filter(testStart:end) - Udc_PSC_filter(testStart:end) ) ./ Udc_PSC_filter(testStart:end) );
Err_Udc = max(delta_Udc)*100

% 三相电感电流abc瞬时值
Iabc_PSC_filter = filter(Num_ac,1,Iabc_PSC);
Iabc_CPP_filter = filter(Num_ac,1,Iabc_CPP);

% 幅值相角计算
for k = 1:size(Iabc_PSC_filter,1)
    [Mag_PSC_filter(k) phase_PSC_filter(k)] = cal_mag_theta(Iabc_PSC_filter(k,:)');
    [Mag_CPP_filter(k) phase_CPP_filter(k)] = cal_mag_theta(Iabc_CPP_filter(k,:)');
end

% delta_Mag = abs( ( Mag_CPP_filter(10000:end) - Mag_PSC_filter(10000:end) ) ./ Mag_PSC_filter(10000:end) );
delta_Mag = abs( ( Mag_CPP_filter(testStart:end) - Mag_PSC_filter(testStart:end) ) ./ Mag_PSC_filter(testStart:end) );
Err_Mag = max(delta_Mag)*100

% 相角差计算
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

% Err_phase = max(abs(delta_phase(10000:end)));
Err_phase = max(abs(delta_phase(testStart:end)))

% 实部虚部计算
Ir_PSC_filter = (2/3)*Iabc_PSC_filter(:,1) - (1/3)*Iabc_PSC_filter(:,2) ...
     - (1/3)*Iabc_PSC_filter(:,3);
Ii_PSC_filter = (1/sqrt(3))*(Iabc_PSC_filter(:,2)-Iabc_PSC_filter(:,3));

Ir_CPP_filter = (2/3)*Iabc_CPP_filter(:,1) - (1/3)*Iabc_CPP_filter(:,2) ...
     - (1/3)*Iabc_CPP_filter(:,3);
Ii_CPP_filter = (1/sqrt(3))*(Iabc_CPP_filter(:,2)-Iabc_CPP_filter(:,3));

TVE = sqrt( ( (Ir_PSC_filter-Ir_CPP_filter).^2+(Ii_PSC_filter-Ii_CPP_filter).^2 ) ./ ( Ir_PSC_filter.^2 + Ii_PSC_filter.^2 ) );
% Err_iabc = max(TVE(10000:end))*100
Err_iabc = max(TVE(testStart:end))*100

%% 作图比较
figure(1)
plot(t,Udc_PSC_filter,t,Udc_CPP_filter);
% 
figure(2)
plot(t,Iabc_PSC_filter,t,Iabc_CPP_filter);
% 
figure(3)
plot(t,Mag_PSC_filter,t,Mag_CPP_filter);
% 
% figure(4)
% plot(t,phase_PSC_filter,t,phase_CPP_filter);
% 
figure(5)
plot(t,delta_phase);
