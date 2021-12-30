% Multimedia Coding
%
%
% AudioStereoToMono.m:  Converts stero audio file to mono 

clear all; close all;

% filenamein = '70stereo.wav';
filenamein = 'DataSet\Audio\70stereo.wav';

info = audioinfo(filenamein);
[x,F] = audioread(filenamein,'native') ; 
fprintf('\n'); 
fprintf('Sampling frequency:      F = %d',F); fprintf(' [Hz] \n'); 
fprintf('Resolution:          nbits = %d',info.BitsPerSample); fprintf(' [bit] \n');

if info.NumChannels == 2
%     xmono = int16(mean(x'));
    xmono = int16(mean(x,2));
end

% filenameout = '70mono.wav';
filenameout = 'DataSet\Audio\70mono.wav';

audiowrite(filenameout,xmono,F);
