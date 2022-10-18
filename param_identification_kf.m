clear all
close all
clc

%JIAYI SU
TODAY = datestr(now)
% This model is from Plett's code, called E2 model
% The tempreture is 25C, parameters with differnet tempreture can be found
% accordingly
% In this model, the noise is coupled with the (input) current
tic

T = 1; % sampling time
soc(1,1) = 0.999; %initial value of soc
ir1(1,1) = 0; % initital value of i_R1
h(1,1) = 0; % initital value of hysteresis
Q = 3.3; % total capacity of the battery cell
C_rate = 1; % C rate (or time used for discharging the cell, 1/2 means uses 2 hours to fully discharge the cell)
num = round((1/C_rate)*3600*soc(1,1)); % # of iterations for the simulation
% I = ones(1,num)*0.4762; % current thru R0, plett
I = ones(1,num)*Q*C_rate; % current thru R0 computed based on C rate, costumized
% ita = 0.9929; % coulombic efficiency
R1 = 0.0025; % R1 in R-C parallel
% C1 = 958.4866; % C in R-C parallel
RC = 4.548; % tau = R1*C1
% RC = 45.48; % tau = R1*C1
R0 = 0.0112;% R0 internal R
gamma = 150; % constant for hysteresis voltage
% ocv(1,1) = 0; % initital value of ocv
K0 = 3.3274; K1=-6.1332e-05; K2=0.0044; K3=0.0498; K4=-0.0107; % ocv regression coeff
M = 0.0443; % hysteresis constant
M0 = 0.0025; %instanies hysteresis constant
s(1,1) = 0;

for k = 1:length(I)
    
    % state of charge x1
    soc(1,k+1) = soc(1,k) - (I(1,k)*T)/(3600*Q) ;
    
    z = soc(1,k);
    
    % current thru r1 x2
    ir1(1,k+1) = exp(-T/(RC))*ir1(1,k) + (1-exp(-T/(RC)))*(I(1,k));
    
    % hysteresis voltage x3
    Ah(1,k) = exp(-abs(((I(1,k))*gamma*T)/(3600*Q)));
    
    h(1,k+1) = Ah(1,k)*h(1,k) + (Ah(1,k) - 1)*sign((I(1,k)));
    
    % open circuit voltage
    ocv(1,k) = K0 - (K1/z) - (K2*z) + (K3*log(z)) + (K4*log(1-z));
    
    % for the S(k) function
    if abs(I(1,k)) > 0
        s(1,k+1) = sign(I(1,k));
    else
        s(1,k+1) = s(1,k);
    end

end

% figure,
% plot(soc)
% 
% figure,
% plot(ir1,'linewidth',2)
% legend('ir1')
% 
% figure,
% plot(h,'linewidth',2)
% legend('hysteresis voltage')
% 
% figure,
% plot(ocv)
% 
% figure,
% plot(soc(1:end-1),ocv)

%%
data = xlsread('data.xlsx',2);
data(1201,:) = [];
data(3453:end,:) = [];
y_real = data(:,end)';

figure,
plot(y_real,'linewidth',2)
hold on
plot(ocv(1:3452),'linewidth',2)
legend('y_k','Estimated OCV')

y_tilde = y_real - ocv(1:3452);
y_tilde_theo = M*h(1:end-1) + M0*s(1:end-1) - R0*I - R1*ir1(1:end-1);

figure,
plot(y_tilde_theo,'linewidth',2)
hold on
plot(y_tilde,'linewidth',2)
legend('ideal y tilde','True y tilde')

%%
A_least = [h(1:length(y_tilde))' s(1:length(y_tilde))' -ir1(1:length(y_tilde))' -I(1:length(y_tilde))'];
% x = lsqnonneg(A,y_tilde')
x_least = A_least\y_tilde'


%% set up kf to estimate the parameters
A_kf = eye(4);
x = NaN(4,length(y_tilde)); %states are unknown parameters 
x(:,1) = [0.01,0.01,0.01,0.01]'; % initial states 
W = 0.001;
w = sqrt(W)*randn(1,length(y_tilde));

xhat(:,1) = [0,0,0,0]';
% xhat(:,1) = [0.1,0.1,0.1,0.1]';

P(:,:,1) = 1*eye(4);

for i = 1:length(y_tilde)
    
    C = A_least(i,:);
    
%     x(:,i+1) = A_kf*x(:,i);
    
%     y(i) = C*x(:,i) + w(i);
    
    K = (P(:,:,i)*C')/(C*P(:,:,i)*C' + W);
    
    xhat(:,i+1) = xhat(:,i) + K*(y_tilde(i) - C*xhat(:,i));
    
    P(:,:,i+1) = P(:,:,i) - (P(:,:,i)*C'*C*P(:,:,i))/(C*P(:,:,i)*C' + W);
end

xx = xhat(:,end)

figure,
plot(1:i,xhat(1,1:end-1),'linewidth',2)
hold on
plot(1:i,xhat(2,1:end-1),'linewidth',2)
hold on
plot(1:i,xhat(3,1:end-1),'linewidth',2)
hold on
plot(1:i,xhat(4,1:end-1),'linewidth',2)
hold on
legend('M','M_0','R_1','R_0')
xlabel('Time index')
xlabel('Variation of parameter values in time')

y_est = xx(1)*h(1:3452) + xx(2)*s(1:3452) - xx(3)*I(1:3452) - xx(4)*ir1(1:3452) + ocv(1:3452);
y_tilde_estimate = y_est - ocv(1:3452);

figure,
plot(y_real,'linewidth',2)
hold on;
plot(y_est,'linewidth',2)
legend('y_k','estimated y_k')

figure,
plot(y_tilde,'linewidth',2)
hold on;
plot(y_tilde_estimate,'linewidth',2)
legend('True y tilde','Estimated y tilde')

MeanSquareError = immse(y_real,y_est)

toc
