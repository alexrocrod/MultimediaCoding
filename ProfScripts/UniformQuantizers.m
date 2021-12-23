% Multimedia Coding
%
%
% UniformQuantizers.m: Examples of midrise and midtread uniform quantizers

clear all; close all;


Xmax = 7.5;

x = -Xmax:.01:Xmax;

Delta = 2;  % Quantization step = 2


% Midrise quantizer
x_mr = floor((double(x)/Delta))*Delta + Delta/2;

figure(1); 
plot([0 0],[-Xmax-1 Xmax+1],'k',[-Xmax-1 Xmax+1],[0 0],'k'); grid on; hold on;
plot(x,x_mr,'r'); axis([-Xmax-1, Xmax+1, -Xmax-1, Xmax+1]); 
title('Midrise Uniform Quantizer'); hold off;



% Midtread quantizer
x_mt = round((double(x)/Delta))*Delta;  

figure(2); 
plot([0 0],[-Xmax-1 Xmax+1],'k',[-Xmax-1 Xmax+1],[0 0],'k'); grid on; hold on;
plot(x,x_mt,'r'); axis([-Xmax-1, Xmax+1, -Xmax-1, Xmax+1]); 
title('Midtread Uniform Quantizer'); hold off;



