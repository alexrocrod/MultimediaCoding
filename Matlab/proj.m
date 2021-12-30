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

if L ==2
    Rs = [2,4];
else
    Rs = [1,2];
end

%% Load Audio
filename = 'Audio\70stereo.wav';

info = audioinfo(filename);
[x,F] = audioread(filename,'native') ; 
fprintf('\n'); 
fprintf('Sampling frequency:      F = %d',F); fprintf(' [Hz] \n'); 
fprintf('Resolution:          nbits = %d',info.BitsPerSample); fprintf(' [bit] \n');

if info.BitsPerSample == 8
    x = int16(x) - 127;
end

%% Convert to mono (needed ??)

if info.NumChannels == 2
    xmono = int16(mean(x,2));
end

%% 











