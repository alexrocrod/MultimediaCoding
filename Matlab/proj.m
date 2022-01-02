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

L = 4; % 2 or 4

if L == 2
    Rs = [2, 4];
elseif L == 4 
    Rs = [1, 2];
else
    disp('Invalid vector dimension L, valid values are 2 or 4.')
    return;
end

Ks = 2.^(Rs.*L);

nbit = 16;



%% Train
% Load Audio
Naud = 64; % 49, 54, 58, 64, 70
[x,F,Nx,maxX] = loadaudio(Naud);

iL = 2;
K = Ks(iL);
R = Rs(iL);
epsilon = 1e-2;

savefile = ['CodeBooks\cb_' num2str(L)  '_' num2str(R) '_' num2str(epsilon) '_' num2str(Naud) '.mat'];
if isfile(savefile)
    disp('Codebook loaded from file.')
    load(savefile)
else
    % Training set 
    N = 1000*K;
    T = zeros(N,L); %% for now
    
    idxT = 1;
    deltai = N*L/2;
    for i = Nx/2-deltai:L:Nx/2+deltai-1
        for j = 1:L
            T(idxT,j) = x(i+j-1); 
        end
        idxT = idxT + 1;
    end
    
    % initial codebook
    tic
    y = split(T,epsilon,K,maxX,L)
    toc    
    save(savefile,'y')
    
    plot(T(:,1),T(:,2),'g.')
    hold on
    plot(y(:,1),y(:,2),'r*')
end
%% Encode
% [x,F,Nx,maxX] = loadaudio(70);
[x,F,Nx,maxX] = loadaudio(1,'music\saynada.wav');
x = x(round(0.25*Nx):round(0.30*Nx));
Nx = length(x);
maxX = max(x);
x1 = zeros(round(Nx/L)+1,L); %% for now

idxT = 1;
for i = 1:L:Nx
    for j = 1:L
        x1(idxT,j) = x(i+j-1); 
    end
    idxT = idxT + 1;
end
tic
[x2,ads] = quantizer(x1,y,L,K);
toc
figure(2)
plot(x1(:,1),x1(:,2),'g.')
hold on
plot(x2(:,1),x2(:,2),'r*')


%% Decode
x_end = int16(zeros(Nx,1)); 

idxT = 1;
for i = 1:L:Nx
    for j = 1:L
         x_end(i+j-1) = x2(idxT,j);
    end
    idxT = idxT + 1;
end

%% Compare 

figure(3)
plot(abs(x_end-x))

sound(single(x)/(2^nbit),F); % original signal
figure(4);
plot(x(152000:154000,:),'g-');
hold on
disp('Press any key to continue');
pause;
clear sound % Added to stop sound playing

sound(single(x_end)/(2^nbit),F); % original signal
plot(x_end(152000:154000,:),'r-');
disp('Press any key to continue');
pause;
clear sound % Added to stop sound playing

%% Functions

function b = split(T,epsilon,K,maxX,L)
    y1 = randi([-maxX maxX],1,L);
    b = y1;
    m = 1;
    while m < K
        n=m;
        m=m+n;
        bnew = [b(1:n,:)+epsilon(ones(n,1),:); b(1:n,:)-epsilon(ones(n,1),:)];
        b = LBG(T,bnew,epsilon);
    end
end

function b = LBG(T,b1,epsilon)
    y = b1;
    [N,~] = size(T);
    [K,~] = size(y);
    D = inf;
    iters = 0;
    while true 
        iters = iters+1;
        % Optimal Partition
        address = zeros(N,1);
        cards = zeros(K,1);
        for i=1:N
            mindist = inf;
            min_idx = 1;
            for j = 1:K
                dist = sqrt((T(i,:)-y(j,:)).^2);
                if dist < mindist
                    mindist = dist;
                    min_idx = j;
                end
            end
            cards(min_idx) = cards(min_idx) + 1;
            address(i) = min_idx;
        end
        
        % New Codebook
        for i=1:K
            if cards(i) == 0 
                [m,cell] = max(cards);
                ys2 = T(address==cell,:);
                y(i,:) = ys2(randi(m));
            else                y(i,:) = sum(T(address==i,:))./cards(i);
            end
        end

        % Evaluate Distortion
        Dold = D;
        D = 0;
        for i=1:N 
            D = D + norm(T(i,:)-y(address(i)))^2;
        end
        D = D/N;
        if abs((Dold-D)/D) < epsilon
            break
        end
            
    end
    b = y;
    fprintf('Ended LBG routine: size=%d, D=%d, iters=%d\n', K, D, iters)
end


function [x2,address] = quantizer(x1,b,L,K)
    N = length(x1(:,1));
    x2 = zeros(N,L);
    address = int16(zeros(N,1));
    cards = int16(zeros(K,1));
    for i=1:N
        i
        mindist = inf;
        min_idx = 1;
        for j = 1:K
            dist = sqrt((x1(i,:)-b(j,:)).^2);
            if dist < mindist
                mindist = dist;
                min_idx = j;
            end
        end
        cards(min_idx) = cards(min_idx) + 1;
        address(i) = min_idx;
        x2(i,:) = b(min_idx,:);
    end
      
end

function [x,F,Nx,maxX] = loadaudio(Naud,file)
    if nargin < 2 
        filename = ['Audio\' num2str(Naud)  'mono.wav'];
    else
        filename = file;
    end
    info = audioinfo(filename);
    [x,F] = audioread(filename,'native') ; 
    fprintf('Sampling frequency:  F = %d [Hz] \n',F); 
    fprintf('Resolution:          nbits = %d [bit] \n',info.BitsPerSample);
    
    % Upscale to 16 bit/sample
    
    if info.BitsPerSample == 8
        disp('Upscaled to 16 bit/sample')
        x = int16(x) - 127;
    end
    
    % Convert to mono
    if info.NumChannels == 2
        disp('Converted to Mono')
        x = int16(mean(x,2));
    end
    Nx = length(x);
    maxX = max(x);
end

