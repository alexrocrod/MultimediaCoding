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

Delta = 2^8;  % quantization with quantization step  2^8
xmono16bits = int16(sign(xmono)).*int16(floor((single(abs(xmono))/Delta))); 
xmono8bits = uint8(xmono16bits + 127);

% filenameout = '70mono8bits.wav';
filenameout = 'DataSet\Audio\70mono8bits.wav';

audiowrite(filenameout,xmono8bits,F,'BitsPerSample',8);