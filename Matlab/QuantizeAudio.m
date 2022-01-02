% Multimedia Coding
%
%
% QuantizeAudio.m:  quantize audio file using a midrise uniform quantizer

clear all; close all;

% filename = '70mono.wav';
filename = 'Audio/70mono.wav';
 
info = audioinfo(filename);
[x,F] = audioread(filename,'native') ;  % load the audio file il native format (integer)
fprintf('\n'); 
fprintf('Sampling frequency:      F = %d',F); fprintf(' [Hz] \n'); 
nbit = info.BitsPerSample;
fprintf('Resolution:          nbits = %d',nbit); fprintf(' [bit] \n \n');

sound(single(x)/(2^nbit),F); % original signal
figure(1); plot(x(152000:154000,:));
disp('Press any key to continue');
pause;
clear sound % Added to stop sound playing

Delta = 2^8;  % quantization with quantization step  2^8
x8 = floor((double(x)/Delta))*Delta + Delta/2;
sound(single(x8)/(2^nbit),F); 
figure(2); plot(x8(152000:154000,:));
disp('Press any key to continue');
pause;
clear sound

Delta = 2^12;  % quantization with quantization step  2^12
x4 = floor((double(x)/Delta))*Delta + Delta/2;
sound(single(x4)/(2^nbit),F); 
figure(3); plot(x4(152000:154000,:));
disp('Press any key to continue');
pause;
clear sound

Delta = 2^13;  % quantization with quantization step  2^13
x3 = floor((double(x)/Delta))*Delta + Delta/2;
sound(single(x3)/(2^nbit),F); 
figure(4); plot(x3(152000:154000,:));
disp('Press any key to continue');
pause;
clear sound

Delta = 2^14;  % quantization with quantization step  2^14
x2 = floor((double(x)/Delta))*Delta + Delta/2;
sound(single(x2)/(2^nbit),F); 
figure(5); plot(x2(152000:154000,:));
disp('Press any key to continue');
pause;
clear sound

Delta = 2^15;  % quantization with quantization step  2^15
x1 = floor((double(x)/Delta))*Delta + Delta/2;
sound(single(x1)/(2^nbit),F); 
figure(6); plot(x1(152000:154000,:));
disp('Press any key to end');
pause;
clear sound

