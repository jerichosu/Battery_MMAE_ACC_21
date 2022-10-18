clear all
close all
clc
%JIAYI SU
TODAY = datestr(now)
fprintf("THIS BATTERY WAS TESTED IN DR.JOSSE'S LAB UNDER THE HELP OF AKILA \n")
fprintf("THE IDENTIFIED PARAMETERS CAN BE FOUND IN identification_coefficients.m \n")
fprintf("THE BATTERY CELL WAS TESTED AT THE ROOM TEMPERATURE (22C) \n")
fprintf("\n")
fprintf("Cell name: LiFePO4 26650 Rechargeable Cell \n")
fprintf("Cell nominal output voltage: 3.2V \n")
fprintf("Cell total capacity: C = 3300mAh (3.3Ah) \n")
fprintf("Other information: Max 10A C rate & 10Wh \n")
tic

T = 1; % sampling time
soc(1,1) = 0.99; %initial value of soc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ir1(1,1) = 0; % initital value of i_R1
h(1,1) = 0; % initital value of hysteresis
Q = 3.3; % total capacity of the battery cell
C_rate = 1/1; % C rate (or time used for discharging the cell, 1/2 means uses 2 hours to fully discharge the cell)
num = round((1/C_rate)*3600*soc(1,1)) - 20; % # of iterations for the simulation
% I = ones(1,num)*0.4762; % current thru R0, plett
I = ones(1,num)*Q*C_rate; % current thru R0 computed based on C rate, costumized
ita = 1; % coulombic efficiency
R1 = 0.0117; % R1 in R-C parallel
% C1 = 958.4866; % C in R-C parallel
RC = 4.548; % tau = R1*C1
R0 = 0.0388;% R0 internal R
gamma = 150; % constant for hysteresis voltage
ocv(1,1) = 0; % initital value of ocv
K0 = 3.3274; K1=-6.1332e-05; K2=0.0044; K3=0.0498; K4=-0.0107; % ocv regression coeff
M = 0.0013; % hysteresis constant
M0 = 0; %instanies hysteresis constant

%initialize the state estimate
zhat = NaN(1,length(I));
hhat = NaN(1,length(I));
ir1hat = NaN(1,length(I));
% zhat(1,1) = 0.8; % INITIAL SOC ESTIAMTE, CHANGE THIS TO SEE THE ESITMATION RESULT!!!
hhat(1,1) = 0;
ir1hat(1,1) = 0;
xhat(:,1) = [zhat(1,1); ir1hat(1,1); hhat(1,1)];
%initialize the output estimate
vbhat = NaN(1,length(I));
s(1,1) = 0;
%**************************************************
%initialize the error covariance matrix
P = NaN(3,3,length(I));
P(:,:,1) = diag([1e-6,1e-8,2e-4]);
% P(:,:,1) = diag([0.01,0.01,0.01]);
% P(:,:,1) = diag([10,10,10]);
%*************************************************

%*************************************************
% initialize the noise
V = 1e-7;
W = 1e-6;
% V = 0.0001;
% W = 0.0001;
%************************************************
v = sqrt(V)*randn(1,length(I));
w = sqrt(W)*randn(1,length(I));

hwait = waitbar(0,'ESTIMATING SOC USING EKF...');
for k = 1:length(I)
    
    % state of charge x1
    soc(1,k+1) = soc(1,k) - ((I(1,k) + v(1,k))*T*ita)/(3600*Q) ;
    
    z = soc(1,k);
    
    % current thru r1 x2
    ir1(1,k+1) = exp(-T/(RC))*ir1(1,k) + (1-exp(-T/(RC)))*(I(1,k) + v(1,k));
    
    % hysteresis voltage x3
    Ah(1,k) = exp(-abs(((I(1,k)+v(1,k))*gamma*T)/3600*Q));
    
    h(1,k+1) = Ah(1,k)*h(1,k) + (Ah(1,k) - 1)*sign((I(1,k) + v(1,k)));
    
    % open circuit voltage
    ocv(1,k) = K0 - (K1/z) - (K2*z) + (K3*log(z)) + (K4*log(1-z));
    
    % for the S(k) function
    if abs((I(1,k) + v(1,k))) > 0
        s(1,k+1) = sign((I(1,k) + v(1,k)));
    else
        s(1,k+1) = s(1,k);
    end
    
    %output voltage y
    vb(1,k) = ocv(1,k) + M*h(1,k) + M0*s(1,k) - (ir1(1,k)*R1) - (I(1,k)*R0) + w(1,k);
    
    
    %%%%%%%%%%%%%%%%%%SIGNAL PROCESSING%%%%%%%%%%%%%%%%%%%%
    %%
    %EKF
    % FOR COMPUTING THE INITIAL SOC ESTIMATE NUMERICALLY
    if k == 1
        
%         %ADDING THE MEASUREMENT NOISE (ASSUME WE CAN ONLY MEASURE THE VOLTAGE)
%         vb1 = ocv(1,k) - ((I(1,k))*R0);
%         fun = @(x) K0 - (K1/x) - (K2*x) + (K3*log(x)) + (K4*log(1-x)) - ((I(1,k) + v(1,k))*R0) - vb1;
%         value = real(fsolve(fun,0.1));
%         %THIS STEP NEEDS TO BE CONSIDERED 
%         zhat(1,k) = min(0.999,value);
%         xhat(1,k) = zhat(1,k)
        zhat(1,k) = 0.5; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        xhat(1,k) = zhat(1,k)
    end
    
    %STEP 1: JACOBIAN MATRIX
    A = [1    0          0;
         0    exp(-T/RC) 0;
         0    0          Ah(1,k)];
    
%     B = [(-T*ita)/(3600*Q)     0;
%         (1 - exp(-T/RC))      0;
%         0                    -abs(gamma*T/(3600*Q))*Ah(1,k)*(1 + sign((I(1,k)+v(1,k)))*hhat(1,k))];
    B = [(-T*ita)/(3600*Q)     0;
        (1 - exp(-T/RC))      0;
        0                    Ah(1,k) - 1];
    
    C = [K1/(zhat(1,k)^2) - K2 + K3/zhat(1,k) + K4/(zhat(1,k)-1) -R1 M];
    
    F = B(:,1);
    G = 1;
    
    %STEP 2:
    % open circuit voltage estimate
    ocvhat(1,k) = K0 - (K1/zhat(1,k)) - (K2*zhat(1,k)) + (K3*log(zhat(1,k))) + (K4*log(1-zhat(1,k)));
    % output voltage estimate
    vbhat(1,k) = ocvhat(1,k) + M*hhat(1,k) + M0*s(1,k) - (ir1hat(1,k)*R1) - (I(1,k)*R0);
    
    % Step 2a: State estimate time update 
    % (FROM PLETT'S CODE, HE LINEARIZED THE MODEL AND USE THIS MODEL TO PRODUCE THE STATE ESTIMATE)
    % xhat(:,k) = A*xhat(:,k) + B*[I; sign(I)];
    
%     WITH NOISE (ASSUME WE CAN ONLY MEASURE THE INPUT SIGNAL)
    f(:,k) = [zhat(1,k) - ((I(1,k) + v(1,k))*T*ita)/(3600*Q);
        exp(-T/(RC))*ir1hat(1,k) + (1-exp(-T/(RC)))*(I(1,k) + v(1,k))
        Ah(1,k)*hhat(1,k) + (Ah(1,k) - 1)*sign((I(1,k) + v(1,k)))];
    
%     %WITHOUT NOISE (ASSUME WE KNOW EXACTLY THE INPUT SIGNAL)
%     f = [zhat(1,k) - (I(1,k)*T*ita)/(3600*Q);
%         exp(-T/(RC))*ir1hat(1,k) + (1-exp(-T/(RC)))*I(1,k); 
%         Ah(1,k)*hhat(1,k) + (Ah(1,k) - 1)*sign(I(1,k))];
    
%         help maintain robustness based on plett, NOT USEFUL FOR HERE (WORSE ESTIMATE) 
%         f(1) = min(1.05,max(-0.05,xhat(1,k))); 
%         f(3) = min(1,max(-1,xhat(3,k)));
    
    %STEP 3: SET UP EKF
    %KALMAN GAIN
    K = (A*P(:,:,k)*C')/(C*P(:,:,k)*C'+G*W*G');
    %STATES ESTIMATE
    xhat(:,k+1) = f(:,k) + K*(vb(1,k) - vbhat(1,k));
    %RICCATI EQUATION
    P(:,:,k+1) = A*P(:,:,k)*A' - K*(C*P(:,:,k)*A') + F*V*F';
    
    % Help maintain robustness based on plett's code
    [~,S,VV] = svd(P(:,:,k+1));
    HH = VV*S*VV';
    P(:,:,k+1) = (P(:,:,k+1) + P(:,:,k+1)' + HH + HH')/4;
    
    %ESTIMATED STATES
    zhat(1,k+1) = xhat(1,k+1);
    ir1hat(1,k+1) = xhat(2,k+1);
    hhat(1,k+1) = xhat(3,k+1);
    
    if mod(k,10)==0,
        waitbar(k/length(I),hwait);
    end;
end
close(hwait);


%% *********************** BANK OF KALMAN FILTER *****************************
% initialize the bank of kalman filters
socbank = 0.01:0.01:0.99; %quantize the soc
numfil = length(socbank); % number of filters
imax = num;
y = vb; % actual output (terminal voltage) 
weight = zeros([numfil,imax+1]); %initialize the weights
cw = zeros([numfil,imax]); %normorlizing factor
weight_soc = 0.99;
weight(1:94,1) = (1-weight_soc)/numfil; % uniform initial guess for the initial weights
weight(96:end,1) = (1-weight_soc)/numfil; % uniform initial guess for the initial weights
weight(95,1) = weight_soc; % uniform initial guess for the initial weights
% yhat = zeros([numfil,imax]); %initialize the estimated output
inno = zeros([numfil,imax]); % initialize the innovation sequence
soc_hat=zeros(1,imax); % initialize the soc estimate using bank of KF
xbankhat = zeros([4,imax+1,numfil]); 
ini_x =[0 0 0 0]';
Pbank=zeros(4,4,imax+1,numfil); %initialize the error covariance matrix
%**********************************************************
% Pt = diag([1e-10 1e-10 1e-10]);
% Pt = diag([1e-2 1e-2 1e-2]);
Pt = diag([1 1 1 1]);
%**********************************************************
soc_hat(1,1)=zhat(1,1); % initial value of bank of kalman filter
%**********************************************************
yyhat=zeros(1,imax);
ocvestbank = zeros(1,imax);
for j=1:numfil
    xbankhat(:,1,j) = ini_x;
    Pbank(:,:,1,j) = Pbank(:,:,1,j) + Pt;
end

pp = 0;

 shat(1,1) = 0;

hwait = waitbar(0,'ESTIMATING SOC USING BANK OF KF...');
for i = 1:imax
    for j = 1:numfil
        
        if pp==1
            weight(:,i) =  1/numfil;
            if j==numfil
                pp=0;
            end
        end
                 
        zz = socbank(j);
%         xbankhat(1,i,j) = zz;
        ocvhatbank(j,i) = K0 - (K1/zz) - (K2*zz) + (K3*log(zz)) + (K4*log(1-zz));
        
        %s term
        if abs((I(1,i) + v(1,i))) > 0
            shat(1,i+1) = sign((I(1,i) + v(1,i)));
        else
            shat(1,i+1) = shat(1,i);
        end
    
        Cbank = [M0 1 -R1 M];
        Dbank = -R0;
%         vbhatbank(j,i) = ocvhatbank(j,i) + M*xbankhat(3,i,j) + M0*s(1,i) - (xbankhat(2,i,j)*R1) - ((I(1,i) + v(1,i))*R0);
        vbhatbank(j,i) = Cbank*[shat(1,i) ocvhatbank(j,i) xbankhat(3,i,j) xbankhat(4,i,j)]' + ((I(1,i) + v(1,i))*Dbank);

%         inno(j,i) = y(i) - Cbank*xbankhat(:,j,i);
        inno(j,i) = y(i) - vbhatbank(j,i);        
        
        Ahbank(1,i) = exp(-abs(((I(1,i)+ v(1,i))*gamma*T)/3600*Q));
            
        Abank = [1     0     0          0;
                 0     1     0          0;
                 0     0     exp(-T/RC) 0;
                 0     0     0          Ahbank(1,i)];
        
        Bbank = [ 0                      0;
                  0                      0;
                  (1 - exp(-T/RC))       0;
                  0                      Ahbank(1,i) - 1];
        Fbank = Bbank(:,1);
         
%         Cbank = [K1/(zz^2) - K2 + K3/zz + K4/(zz-1) R1 M];
        
        %KALMAN GAIN
        Kbank = (Abank*Pbank(:,:,i,j)*Cbank')/(Cbank*Pbank(:,:,i,j)*Cbank'+G*W*G');
        %STATES ESTIMATE
        xbankhat(:,i+1,j) = Abank*xbankhat(:,i,j) + Bbank*[I(1,i); sign(I(1,i))] + Kbank*inno(j,i);
        %RICCATI EQUATION
        Pbank(:,:,i,j) = (Abank - Kbank*Cbank)*Pbank(:,:,i,j)*(Abank - Kbank*Cbank)' - Kbank*(Cbank*Pbank(:,:,i,j)*(Abank - Kbank*Cbank)') + Fbank*V*Fbank';
        
          % Weight Update Equations
        S = Cbank*Pbank(:,:,i+1,j)*Cbank' + G*W*G';

        % Weight before normalization
        cw(j,i) = ((abs(S))^(-0.5))*(exp(((-0.5)*(inno(j,i))*(inno(j,i)))/(S)))*weight(j,i);
        %cw = ((2*pi)^(-0.5))*((abs(S))^(-0.5))*(exp(((-0.5)*(inno)*(inno))/(S)))*w;
        %[x(:,i+1,j),P(:,:,i+1,j),yhat(j,i),cw(j,i),inno(j,i)] = soc_kalman_filter(x(:,i,j),P(:,:,i,j),y(i),I,SOC_est(j),W,V,w(j,i));    
        
    end
    c = sum(cw(:,i));
    weight(:,i+1) = cw(:,i)./c;
    
    [Y, II] = max(weight(:,i+1));
    soc_hat(1,i+1) = socbank(II);
    yyhat(1,i) = vbhatbank(II,i);
    ocvestbank(1,i) = ocvhatbank(II,i);
    pp=pp+1;
    
    if mod(i,10)==0,
        waitbar(i/length(I),hwait);
    end;
end
close(hwait);

SOC_AE_BANK = abs(soc - soc_hat)*100;
SOC_AE_EKF = abs(soc - zhat)*100;

% figure,
% plot(1:length(I)+1,soc_hat,1:length(I)+1, soc,'linewidth', 2)
% legend('SOC estimate using BANK OF KF','SOC')
% title('SOC')
% xlabel('Time (seconds), k')
% ylabel('STATE OF CHARGE ESTIMATE')

figure,
set(gcf, 'position', [200 300 600 650]);
subplot(2,1,1);
plot( 1:length(I)+1, zhat,'r', 1:length(I)+1, soc_hat,'b',1:length(I)+1, soc,'y', 'linewidth', 2.2)
legend('SOC Estimate Using EKF','SOC Estimate Using MMAE','True SOC')
title('(a) SOC Estimation Results at 1C Rate')
xlabel('Time k, seconds')
ylabel('State of Charge')
grid on

subplot(2,1,2);
plot(1:length(I)+1, SOC_AE_BANK,'b',1:length(I)+1, SOC_AE_EKF,'r', 'linewidth', 2.2)
legend('Absolute Error Using MMAE','Absolute Error Using EKF')
title('(b) Absolute Error, %')
xlabel('Time k, seconds')
ylabel('Absolute Error, %')
grid on

% figure,
% plot(1:length(I)+1, ir1, 1:length(I)+1, ir1hat,'--','linewidth', 2)
% legend('Current thru RC parallel','Estimated I-R1')
% title('Current thru RC parallel')
% xlabel('Time (seconds), k')
% ylabel('Current thru RC parallel')
% 
% figure,
% plot(1:length(I)+1, h, 1:length(I)+1, hhat,'--','linewidth', 2)
% legend('Hysteresis voltage','Hysteresis voltage estimate')
% title('Hysteresis voltage')
% xlabel('Time (seconds), k')
% ylabel('HYSTERESIS VOLTAGE')

figure,
plot(1:length(I), ocv, 1:length(I), ocvhat,'--',1:length(I),ocvestbank,'linewidth', 2)
legend('Open Circuit Voltage OCV','OCV Estimate using EKF','OCV Estimate using bank of KF')
title('Open Circuit Voltage OCV')
xlabel('Time (seconds), k')
ylabel('OPEN CIRCUIT VOLTAGE')
grid on

figure,
plot(1:length(I), vb, 1:length(I), yyhat,'--', 1:length(I), vbhat,'linewidth', 2)
legend('Output Voltage Vb','Vb Estimate using Bank of KF','Vb Estimate using EKF')
title('Output Voltage Vb')
xlabel('Time (seconds), k')
ylabel('OUTPUT VOLTAGE')
grid on

soc_mse_ekf = immse(soc,zhat)
soc_mse_bkf = immse(soc,soc_hat)

% figure,
% plot(1:length(I),(I + v), 'linewidth', 2)
% legend('Load current I with noise')
% title('Load current I with noise')
% xlabel('Time (seconds), k')
% ylabel('Load current I with noise')

toc



