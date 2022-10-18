clear all
close all
clc

data = xlsread('rc_constant_room_temp.xlsx',1);
A = [data(:,1).*60 data(:,14)];
A(1201,:) = [];
x = A(:,1);
y = A(:,2);
figure,
plot(x,y,'linewidth',2); hold on;
xlim([0 length(x)])
ylim([3.1 3.34])
legend('Raw data for RC constant')
yline(3.3181,'--',{'STEADY STATE'});

five_rc = min(find(y == 3.3181)) - 20*60;
tau = five_rc/5

