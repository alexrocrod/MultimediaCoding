function [x,F,Nx,maxX] = loadaudio(Naud,file,L,K)
    if Naud ~= 1 
        filename = ['Audio\' num2str(Naud)  'mono.wav'];
    else
        filename = file;
    end
    info = audioinfo(filename);
    [x,F] = audioread(filename,'native') ; 
%     fprintf('Sampling frequency:  F = %d [Hz] \n',F); 
%     fprintf('Resolution:          nbits = %d [bit] \n',info.BitsPerSample);
        
    % Convert to mono
    if info.NumChannels == 2
%         disp('Converted to Mono')
        x = int16(mean(x,2));
    end
    Nx = length(x);
    if nargin > 2 && round(Nx/L)-1 > 2000 * K 
%         disp('reducing size')
        delta = 1000 * K;
        mid = round(Nx/2);
        x = x(mid-delta:mid+delta);
        Nx = length(x);    
    end
    maxX = max(x);
end