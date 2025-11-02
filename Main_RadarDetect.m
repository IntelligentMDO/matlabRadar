
clear all;
clc;
close all;

Num_Parallel = 16;                                                           %并行线程数

p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= Num_Parallel
    delete(p);
    parpool(Num_Parallel);
end


%% 加载数据
load('PostData_zFlyingWing.mat')
Data_Post = PostData_zFlyingWing;


%% 定义参数
Act_Solver = [1,1];                                                         %是否激活求解器([前向,后向])(否_0,是_1)

%角域
Mode_DomainPDF = 2;                                                         %角域PDF模式(常量_1,高斯分布_2)
Theta_Total_Divide = {90};                                                  %中心角_Theta
Psi_Total_Divide = {0:5:180};                                               %中心角_Psi
Delta_Theta = 1;                                                            %角域扰动_Theta
Delta_Psi = 3;                                                              %角域扰动_Psi
DReSample = 0.001;                                                          %角度重采样间隔
K_BW = 0.1;                                                                 %核密度估计带宽平滑系数

%雷达
Mode_RadarEq = 2;                                                           %雷达方程模式(定距离_1,定信噪比_2)(影响前向)
if Mode_RadarEq == 1                                                        %***定距离***
    R = 50e3;                                                               %距离(m)
    SNR_Fix = NaN; 

elseif Mode_RadarEq == 2                                                    %***定信噪比***
    R = NaN;
    SNR_Fix = 10*log10(10^(13/10)*2);                                       %峰值信噪比(dB)
end

Mode_RadarPara = 2;                                                         %雷达参数模式(指定_1,归一化_2)(影响前向&后向)
if Mode_RadarPara == 1                                                      %***指定***
    P_t = 1.5e6;                                                            %峰值功率(W)
    G = 45;                                                                 %天线增益(dB)
    Freq = 5.6e9;                                                           %频率(Hz)
    B = 1/0.2e-6;                                                           %带宽(Hz)
    F = 3;                                                                  %噪声系数(dB)
    L = 6;                                                                  %损耗(dB)
elseif Mode_RadarPara == 2
    P_t = NaN;
    G = NaN;
    Freq = NaN;
    B = NaN;
    F = NaN;
    L = NaN;
end

P_fa = 1e-6;                                                                %虚警率
if Act_Solver(1) == 1
    N_Total_Forward = {[0,1]};                                              %脉冲数({[相关,不相关],[相关,不相关],...})
end
if Act_Solver(2) == 1
    N_Total_Backward = {[10,0],[0,10]};                                     %脉冲数({[相关,不相关],[相关,不相关],...})
    P_d_Fix = 0.9;                                                          %给定检测概率
end


%% 预处理
fprintf('\nRun Pretreat...');

Theta_Total = [Theta_Total_Divide{:}];
Psi_Total = [Psi_Total_Divide{:}];

global RunNum_Local RunNum_Total;
RunNum_Local = 0;
RunNum_Total = length(Theta_Total) * length(Psi_Total);
Q = parallel.pool.DataQueue;
afterEach(Q, @UpdateParforProgress);

% tic
fprintf(' ');
for nTheta = 1:length(Theta_Total)
    parfor nPsi = 1:length(Psi_Total)

        Theta_Center = Theta_Total(nTheta);
        Psi_Center = Psi_Total(nPsi);

        [Mean_RCS_m2_Domain_Total(nTheta,nPsi),...
         R_Total(nTheta,nPsi),...
         SNR_Total(nTheta,nPsi),...
         Gamma_dBsm_KDE_Total{nTheta,nPsi},...
         p_Gamma_dBsm_KDE_Total{nTheta,nPsi}]...
        = RadarDetect_GenPDF(Data_Post, DReSample, K_BW,...
                             Mode_DomainPDF, Theta_Center, Psi_Center, Delta_Theta, Delta_Psi,...
                             Mode_RadarEq, R, SNR_Fix,...
                             Mode_RadarPara, P_t, G, Freq, B, F, L);
        
        send(Q, (nTheta-1)*length(Psi_Total)+nPsi);
        
    end
end

for nTheta = 1:length(Theta_Total)
    parfor nPsi = 1:length(Psi_Total)
        
        SNR_Local = SNR_Total(nTheta, nPsi);
        Gamma_dBsm_KDE_Local = Gamma_dBsm_KDE_Total{nTheta, nPsi};
        p_Gamma_dBsm_KDE_Local = p_Gamma_dBsm_KDE_Total{nTheta, nPsi};
        
        F_p_Gamma_Forward_Total{nTheta,nPsi}...
        = @(Gamma) interp1(Gamma_dBsm_KDE_Local, p_Gamma_dBsm_KDE_Local, Gamma, 'pchip',0);

        IntRange_Gamma_Forward_Total{nTheta,nPsi}...
        = Gamma_dBsm_KDE_Local([1,end]);

        F_p_Gamma_Backward_Total{nTheta,nPsi}...
        = @(Gamma, SNR) interp1(Gamma_dBsm_KDE_Local - SNR_Local + SNR, p_Gamma_dBsm_KDE_Local, Gamma, 'pchip',0);

        IntRange_Gamma_Backward_Total{nTheta,nPsi}...
        = @(SNR) Gamma_dBsm_KDE_Local([1,end]) - SNR_Local + SNR
        
    end
end

fprintf('\n');
% toc


%% 前向模式_单脉冲检测概率
% tic
if Act_Solver(1) == 1
    fprintf('Run Solver_Forward...');

    for nN = 1:length(N_Total_Forward)
        
        global RunNum_Local RunNum_Total;
        RunNum_Local = 0;
        RunNum_Total = length(Theta_Total) * length(Psi_Total);
        Q = parallel.pool.DataQueue;
        afterEach(Q, @UpdateParforProgress);

        clear P_d_Total_Cache;
            
        fprintf(' ');
        for nTheta = 1:length(Theta_Total)
            parfor nPsi = 1:length(Psi_Total)
                
                F_p_Gamma_Local = F_p_Gamma_Forward_Total{nTheta, nPsi};
                IntRange_Gamma_Local = IntRange_Gamma_Forward_Total{nTheta, nPsi};
                N_Total_Local = N_Total_Forward{nN};

                P_d_Total_Cache(nTheta,nPsi) = RadarDetect_Int(F_p_Gamma_Local, IntRange_Gamma_Local, P_fa, N_Total_Local);

                send(Q, (nTheta-1)*length(Psi_Total)+nPsi);
                
            end
        end

        P_d_Total{nN} = P_d_Total_Cache;
        
    end
    
    fprintf('\n');
end
% toc


%% 后向模式_多脉冲积累增益
if Act_Solver(2) == 1
    fprintf('Run Solver_Backward...');
    
    % tic
    %单脉冲
    global RunNum_Local RunNum_Total;
    RunNum_Local = 0;
    RunNum_Total = length(Theta_Total) * length(Psi_Total);
    Q = parallel.pool.DataQueue;
    afterEach(Q, @UpdateParforProgress);

    fprintf(' ');
    for nTheta = 1:length(Theta_Total)
        parfor nPsi = 1:length(Psi_Total)

            F_p_Gamma_Local = F_p_Gamma_Backward_Total{nTheta, nPsi};
            IntRange_Gamma_Local = IntRange_Gamma_Backward_Total{nTheta, nPsi};
            N_Total_Local = 1;

            F_Err = @(SNR) RadarDetect_Int(@(Gamma)F_p_Gamma_Local(Gamma, SNR), IntRange_Gamma_Local(SNR), P_fa, N_Total_Local)...
                           - P_d_Fix;
            
            SNR_Lb = 10*log10(10^(13/10)*2) - 10;
            SNR_Ub = 10*log10(10^(13/10)*2) + 20;

            [SNR_Require_SinglePulse_Total(nTheta,nPsi), Error_SinglePulse_Total(nTheta,nPsi)]...
            = RadarDetect_BinarySearch(F_Err, SNR_Lb, SNR_Ub, P_d_Fix);

            send(Q, (nTheta-1)*length(Psi_Total)+nPsi);
            
        end
    end
    % toc
    
    % tic
    %多脉冲
    global RunNum_Local RunNum_Total;
    RunNum_Local = 0;
    RunNum_Total = length(N_Total_Backward) * length(Theta_Total) * length(Psi_Total);
    Q = parallel.pool.DataQueue;
    afterEach(Q, @UpdateParforProgress);
    
    fprintf(' ');
    parfor nRun = 1:RunNum_Total

        nN = floor((nRun-1)/length(Theta_Total)/length(Psi_Total)) + 1;
        nTheta = floor(mod((nRun-1)/length(Psi_Total), length(Theta_Total))) + 1;
        nPsi = mod(nRun-1, length(Psi_Total)) + 1;
        
        F_p_Gamma_Local = F_p_Gamma_Backward_Total{nTheta, nPsi};
        IntRange_Gamma_Local = IntRange_Gamma_Backward_Total{nTheta, nPsi};
        N_Total_Local = N_Total_Backward{nN};

        F_Err = @(SNR) RadarDetect_Int(@(Gamma)F_p_Gamma_Local(Gamma, SNR), IntRange_Gamma_Local(SNR), P_fa, N_Total_Local)...
                       - P_d_Fix;

        SNR_Lb = 10*log10(10^(SNR_Require_SinglePulse_Total(nTheta, nPsi)/10) / sum(N_Total_Local) / 2);
        SNR_Ub = SNR_Require_SinglePulse_Total(nTheta, nPsi);
        
        [SNR_Require_MultiPulse_Total_Cache_1(nRun), Error_MultiPulse_Total_Cache_1(nRun)]...
        = RadarDetect_BinarySearch(F_Err, SNR_Lb, SNR_Ub, P_d_Fix);

        send(Q, nRun);
    end

    for nRun = 1:RunNum_Total  

        nN = floor((nRun-1)/length(Theta_Total)/length(Psi_Total)) + 1;
        nTheta = floor(mod((nRun-1)/length(Psi_Total), length(Theta_Total))) + 1;
        nPsi = mod(nRun-1, length(Psi_Total)) + 1;

        SNR_Require_MultiPulse_Total_Cache_2(nN,nTheta,nPsi) = SNR_Require_MultiPulse_Total_Cache_1(nRun);
        
        Error_MultiPulse_Total_Cache_2(nN,nTheta,nPsi) = Error_MultiPulse_Total_Cache_1(nRun);

    end

    for nN = 1:length(N_Total_Backward)

        SNR_Require_MultiPulse_Total{nN}...
        = reshape(SNR_Require_MultiPulse_Total_Cache_2(nN,:,:), length(Theta_Total), length(Psi_Total));
        
        Error_MultiPulse_Total{nN}...
        = reshape(Error_MultiPulse_Total_Cache_2(nN,:,:), length(Theta_Total), length(Psi_Total));

    end
    % toc

    %计算距离及增益
    R_SinglePulse_Total = ((R_Total.^4) .*  10.^(SNR_Total/10) ./ 10.^(SNR_Require_SinglePulse_Total/10)) .^ (1/4);

    for nN = 1:length(N_Total_Backward)
        
        R_MultiPulse_Total{nN} = ((R_Total.^4) .*  10.^(SNR_Total/10) ./ 10.^(SNR_Require_MultiPulse_Total{nN}/10)) .^ (1/4);
        
        IntGain_Total{nN} = 10.^(SNR_Require_SinglePulse_Total/10) ./ 10.^(SNR_Require_MultiPulse_Total{nN}/10);

    end
    
    fprintf('\n');
end


%% 赋值并保存
%赋值
Data_RadarDetect.Data_Post = Data_Post;

Data_RadarDetect.Act_Solver = Act_Solver;

Data_RadarDetect.Domain.Mode_DomainPDF = Mode_DomainPDF;
Data_RadarDetect.Domain.Theta_Divide = Theta_Total_Divide;
Data_RadarDetect.Domain.Psi_Divide = Psi_Total_Divide;
Data_RadarDetect.Domain.Theta = Theta_Total;
Data_RadarDetect.Domain.Psi = Psi_Total;
Data_RadarDetect.Domain.Delta_Theta = Delta_Theta;
Data_RadarDetect.Domain.Delta_Psi = Delta_Psi;
Data_RadarDetect.Domain.DReSample = DReSample;
Data_RadarDetect.Domain.K_BW = K_BW;

Data_RadarDetect.RadarEq.Mode_RadarEq = Mode_RadarEq;
Data_RadarDetect.RadarEq.R = R;
Data_RadarDetect.RadarEq.SNR_Fix = SNR_Fix;

Data_RadarDetect.RadarPara.Mode_RadarPara = Mode_RadarPara;
Data_RadarDetect.RadarPara.P_t = P_t;
Data_RadarDetect.RadarPara.G = G;
Data_RadarDetect.RadarPara.Freq = Freq;
Data_RadarDetect.RadarPara.B = B;
Data_RadarDetect.RadarPara.F = F;
Data_RadarDetect.RadarPara.L = L;

Data_RadarDetect.P_fa = P_fa;

Data_RadarDetect.Result_Pre.Mean_RCS_m2_Domain = Mean_RCS_m2_Domain_Total;
Data_RadarDetect.Result_Pre.R = R_Total;
Data_RadarDetect.Result_Pre.SNR = SNR_Total;
Data_RadarDetect.Result_Pre.Gamma_dBsm_KDE_Total = Gamma_dBsm_KDE_Total;
Data_RadarDetect.Result_Pre.p_Gamma_dBsm_KDE_Total = p_Gamma_dBsm_KDE_Total;

if Act_Solver(1) == 1
    for nN = 1:length(N_Total_Forward)
        Data_RadarDetect.Result_Forward(nN).N_Total_Forward = N_Total_Forward{nN};
        Data_RadarDetect.Result_Forward(nN).R = R_Total;
        Data_RadarDetect.Result_Forward(nN).SNR = SNR_Total;
        Data_RadarDetect.Result_Forward(nN).P_d = P_d_Total{nN};
    end
end

if Act_Solver(2) == 1
    for nN = 1:length(N_Total_Backward)
        Data_RadarDetect.Result_Backward(nN).N_Total_Backward = N_Total_Backward{nN};
        Data_RadarDetect.Result_Backward(nN).P_d_Fix = P_d_Fix;
        Data_RadarDetect.Result_Backward(nN).R_SinglePulse = R_SinglePulse_Total;
        Data_RadarDetect.Result_Backward(nN).R_MultiPulse = R_MultiPulse_Total{nN};
        Data_RadarDetect.Result_Backward(nN).SNR_Require_SinglePulse = SNR_Require_SinglePulse_Total;
        Data_RadarDetect.Result_Backward(nN).SNR_Require_MultiPulse = SNR_Require_MultiPulse_Total{nN};
        Data_RadarDetect.Result_Backward(nN).IntGain = IntGain_Total{nN};
        Data_RadarDetect.Result_Backward(nN).Error_SinglePulse = Error_SinglePulse_Total;
        Data_RadarDetect.Result_Backward(nN).Error_MultiPulse = Error_MultiPulse_Total{nN};
    end
end

% %保存
% FileName_Save_Mat = 'Data_RadarDetect';
% eval([FileName_Save_Mat, '=Data_RadarDetect;']);
% save([pwd,'\',FileName_Save_Mat,'.mat'], FileName_Save_Mat);


%% 作图
if Act_Solver(1) == 1

    subplot(5,1,1);
    plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Forward(1).P_d);
    ylabel('P_d');
    ylim([0,1]);
    hold on;
    grid on;
    
    subplot(5,1,2);
    if Mode_RadarEq == 1
        plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Forward(1).SNR);
        ylabel('SNR');
    elseif Mode_RadarEq == 2
        plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Forward(1).R);
        ylabel('R');
    end
    hold on;
    grid on;

end

if Act_Solver(2) == 1
    
    subplot(5,1,3);
    plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Backward(1).IntGain);
    hold on;
    plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Backward(2).IntGain);
    ylabel('IntGain');
    grid on;

    subplot(5,1,4);
    plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Backward(1).R_SinglePulse);
    ylabel('R\_SinglePulse');
    hold on;
    grid on;

    subplot(5,1,5);
    plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Backward(1).R_MultiPulse);
    hold on;
    plot(Data_RadarDetect.Domain.Psi, Data_RadarDetect.Result_Backward(2).R_MultiPulse);
    ylabel('R\_MultiPulse');
    grid on;

end


%%
function UpdateParforProgress(~)

    global RunNum_Local RunNum_Total;

    Num_Str = length(num2str(RunNum_Total));
    if RunNum_Local ~= 0
        fprintf([repmat('\b',1,Num_Str),'\b',repmat('\b',1,Num_Str)]);
    end
    RunNum_Local = RunNum_Local + 1;
    
    fprintf(['%',num2str(Num_Str),'d/%',num2str(Num_Str),'d'], RunNum_Local, RunNum_Total);

end







































