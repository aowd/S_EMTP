function [ A, theta ] = cal_mag_theta( xabc )
%% 本函数用于计算三相对称信号的幅值和相角
% xabc: 输入三相对称信号，abc坐标下的瞬时值，列向量
% A: 输入信号的幅值
% theta: 输入信号的相角

%% 将信号变换到alpha_beta坐标系下
% alpha_beta变换矩阵
Tabc2ab = [
    2/3 -1/3 -1/3
    0 1/sqrt(3) -1/sqrt(3)
    ];

% ab变换： abc ==> alpha_beta
xab = Tabc2ab * xabc;

%% 计算幅值和相角
% 幅值
A = sqrt( xab(1)^2 + xab(2)^2 );

% 相角
if abs(xab(1)) < 1e-10 && abs(xab(2)) < 1e-10
    theta = 0;
elseif abs(xab(1)) < 1e-10
    if xab(2) > 0
        theta = pi/2;
    else
        theta = -pi/2;
    end
elseif abs(xab(2)) < 1e-10
    if xab(1) > 0
        theta = 0;
    else
        theta = pi;
    end
else
    theta = atan( xab(2) / xab(1) );
    if xab(1) < 0
        if theta < 0
            theta = theta + pi;
        else
            theta = theta - pi;
        end
    end
end

end

