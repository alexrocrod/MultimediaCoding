%% Project 3 -  Multimedia Coding
% Name: Alexandre da Rocha Rodrigues
% Matricola: 2039952
% January 2022

close all
clear all

%% Project 3
% CD-quality audio 
% 16bit/sample
% Vector Quantizier - LBG-split algorithm
% Vectors |    L=2    |    L=4    |
% Rates   | R = [2,4] | R = [1,2] | (bit/sample)
% Various samples of voice and music

%% Parameters

L = 2; % 2 or 4

if L == 2
    Rs = [2, 4];
elseif L == 4 
    Rs = [1, 2];
else
    disp('Invalid vector dimension L, valid values are 2 or 4.')
    return;
end

%% Load Audio
filename = 'Audio\70stereo.wav';

info = audioinfo(filename);
[x,F] = audioread(filename,'native') ; 
fprintf('Sampling frequency:  F = %d [Hz] \n',F); 
fprintf('Resolution:          nbits = %d [bit] \n',info.BitsPerSample);

%% Upscale to 16 bit/sample

if info.BitsPerSample == 8
    disp('Upscaled to 16 bit/sample')
    x = int16(x) - 127;
end

%% Convert to mono (needed ??)

if info.NumChannels == 2
    disp('Converted to Mono')
    xmono = int16(mean(x,2));
end

%% Code








%% Compare 








%% Decode















