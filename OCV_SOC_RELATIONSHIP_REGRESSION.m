clear all
close all
clc

ocv_soc_data = xlsread('OCVSOC_data.xlsx',1);
a = ocv_soc_data(:,[1 3 7]);
a(:,1) = a(:,1).*60;
time = a(:,1);
ocvcharge = a(:,2);
ocvdischarge = a(:,3);
ocvmean = (ocvcharge + ocvdischarge)/2;

%ocv vs time
figure,
plot(time,ocvcharge,'linewidth', 2)
hold on
plot(time,ocvdischarge,'linewidth', 2)
hold on
plot(time,ocvmean,'--','linewidth', 2)
legend('True OCV when charging','True (FLIPPED) OCV when discharging','Mean of true OCV')
xlabel('TIME(SECONDS)-- TESTED FOR 200mins')
ylabel('OCV(VOLTAGE)')
title('OCV VS TIME (22 celsius)')
% ocv vs soc
soc = a(:,1)./12000;
figure,
plot(soc,ocvcharge,'linewidth', 2)
hold on
plot(soc,ocvdischarge,'linewidth', 2)
hold on
plot(soc,ocvmean,'--','linewidth', 2)
legend('True OCV when charging','True OCV when discharging','Mean of true OCV')
xlabel('SOC')
ylabel('OCV(VOLTAGE)')
title('OCV VS SOC (22 celsius)')
%%
% param identification
A_param = NaN(length(soc),5);
b = ocvmean;

for j = 1
    for i = 1:length(soc)
        A_param(i,j) = 1;
    end 
end

for j = 2
    for i = 1:length(soc)
        A_param(i,j) = -1/soc(i);
    end 
end

for j = 3
    for i = 1:length(soc)
        A_param(i,j) = -soc(i);
    end 
end

for j = 4
    for i = 1:length(soc)
        A_param(i,j) = log(soc(i));
    end 
end

for j = 5
    for i = 1:length(soc)
        A_param(i,j) = log(1-soc(i));
    end 
end

%linear regression
x = A_param\b;

K0 = x(1)
K1 = x(2)
K2 = x(3)
K3 = x(4)
K4 = x(5)

%% compare to the test data

t = length(soc);
T=1/360;
Q = 3.3;
SOC(1,1)=1;
C_rate = 1/(200/60); % C rate 
num = round((1/C_rate)*inv(T)*SOC(1,1)); % # of iterations for the simulation
I = ones(1,num)*Q*C_rate; % current thru R0 computed based on C rate, costumized
k=length(I);

for i=1:k
    
    SOC(1,i+1)=SOC(1,i)-((I(1,i)*T)/Q);
    z = SOC(1,i);
    % Open Circuit Voltage
    OCV(1,i) = K0 - (K1/z) - (K2*z) + (K3*log(z)) + (K4*log(1-z));
        
end
SOC = SOC(:,2:end-1);
OCV = OCV(:,2:end);
figure,
plot(soc(2:length(SOC(:,2:end-1))+1,:),ocvmean(2:length(SOC(:,2:end-1))+1,:),'b','linewidth', 2.4)
hold on
plot(SOC,OCV,'--r','linewidth', 2.4)
xlabel('State of Charge')
ylabel('Open Circuit Voltage, V')
legend('Measured OCV (Raw data)','Estimated OCV')
title('OCV-SOC Relationship at 22 Celsius')