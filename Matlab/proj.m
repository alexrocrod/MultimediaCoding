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

Ks = 2.^(Rs.*L)

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
idx = 1;
K = Ks(idx);
R = Rs(idx);

epsilon = 0.001;

% Training set ??
T = rand(N,K); %% for now

% Initial Codebook
y = split(T,epsilon);







%% Compare 



% sound(single(x4)/(2^nbit),F); 
% figure(3); plot(x4(152000:154000,:));
% disp('Press any key to continue');
% pause;
% clear sound


%% Decode







%% Functions

function b = split(T,epsilon)
    N,K = size(T);
    y1 = mean(T);
    b = y1;
    while length(b) < K
        bnew = zeros(length(b)*2,N);
        for i=1:length(b)
            bnew = [ bnew ; (1-epsilon)*b(i) ; (1+epsilon)*b(i)];
        end
        b = LBG(bnew,epsilon);
    end
end

function b = LBG(b1,epsilon)
    D = inf;
    Dold = inf;
    
    y = b1;
    while (Dold-D)/D < epsilon

        % Optimal Partition
        for i=1:K
            part(i) = 0; %??
        end
        
        % New Codebook
        for i=1:K
            y(i) = sum(part(i))/length(part(i));
        end
        
        % Evaluate Distortion
        Dold = D;
        D = 0;
        for i=1:N 
            D = D + norm(x-Q(x))^2;
        end
        D = D/N;
            
    end
    b = y;
end













