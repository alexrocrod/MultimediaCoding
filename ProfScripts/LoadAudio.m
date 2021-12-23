% Multimedia Coding
%
%
% LoadAudio.m:  load audio file 

clear all; close all;

% filename = '70stereo.wav';
filename = 'DataSet\Audio\70stereo.wav';

info = audioinfo(filename);
[x,F] = audioread(filename,'native') ; 
fprintf('\n'); 
fprintf('Sampling frequency:      F = %d',F); fprintf(' [Hz] \n'); 
fprintf('Resolution:          nbits = %d',info.BitsPerSample); fprintf(' [bit] \n');

if info.BitsPerSample == 8
    x = int16(x) - 127;
end
