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
epsilon = 1e-2;
L = 2; % 2 or 4
iL = 1; % index to identify the rate selected

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

%% Load Training Audio
% Load Audio
% 1 audio file
% Naud = 64; % 49, 54, 58, 64, 70
% [x,F,Nx,maxX] = loadaudio(Naud);

% All files together
% [x,F,Nx,maxX] = loadAllAudio(); Naud = 100;

% music
% [x,F,Nx,maxX] = loadaudio(1,'music\SayNada.wav'); Naud = 200;

% All audio and music
[x,F,Nx,maxX] = loadAllAudioMusic(); Naud = 300;

% x = x(round(0.25*Nx):round(0.30*Nx));
% music2
% Naud = 210;
% [x,F,Nx,maxX] = loadaudio(1,'music\Good4U.wav');
% x = x(round(0.25*Nx):round(0.30*Nx));
% music3
% Naud = 220;
% [x,F,Nx,maxX] = loadaudio(1,'music\WaitingOnAWar.wav');

%% Train
K = Ks(iL);
R = Rs(iL);

savefile = ['CodeBooks\cb_' num2str(L)  '_' num2str(R) '_' num2str(epsilon) '_' num2str(Naud) '.mat'];
if isfile(savefile)
    fprintf('Codebook loaded from file %s\n',savefile)
    load(savefile)
else
    % Training set 
    N = 1000*K;
    T = zeros(N,L); %% for now
    
    idxT = 1;
    deltai = N*L/2;
    for i = round(Nx/2-deltai):L:round(Nx/2+deltai-1)
        for j = 1:L
            T(idxT,j) = x(i+j-1); 
        end
        idxT = idxT + 1;
    end
    
    % initial codebook
    tic
    y = split(T,epsilon,K,maxX,L);
    toc   
%     fprintf('Codebook: \n')
%     disp(y)
    save(savefile,'y')
    
    plot(T(:,1),T(:,2),'g.')
    hold on
    plot(y(:,1),y(:,2),'r*')
    legend('Original','Codebook')
end

%% Load Encoding Audio
% audio file 
% NEnc = 70; [x,F,Nx,maxX] = loadaudio(NEnc);

% music
% [x,F,Nx,~] = loadaudio(1,'music\SayNada.wav'); NEnc = 200;
[x,F,Nx,~] = loadaudio(1,'music\Good4U.wav'); NEnc = 210;
% [x,F,Nx,~] = loadaudio(1,'music\WaitingOnAWar.wav'); NEnc = 220;
% [x,F,Nx,~] = loadaudio(1,'music\ForeverAfterAll.wav'); NEnc = 230;
% [x,F,Nx,~] = loadaudio(1,'music\SummerThing.wav'); NEnc = 240;
% [x,F,Nx,~] = loadaudio(1,'music\TodoDeTi.wav'); NEnc = 250;
% [x,F,Nx,~] = loadaudio(1,'music\TooOfficial.wav'); NEnc = 260;

% x = x(round(0.25*Nx):round(0.30*Nx));
% Nx = length(x);
% maxX = max(x);

savepic = ['Results\cb_' num2str(L)  '_' num2str(R) '_' num2str(epsilon) '_' num2str(Naud) '\Enc' num2str(NEnc) '\'];
mkdir(savepic(1:end-1))
f = gcf;
exportgraphics(f,[savepic '1.png'])

%% Encode
tic
x1 = zeros(round(Nx/L),L);
idxT = 1;
for i = 1:L:Nx-L
    for j = 1:L
        x1(idxT,j) = x(i+j-1); 
    end
    idxT = idxT + 1;
end

[x2,ads,D] = quantizer(x1,y,L,K);
toc

fprintf('Audio Encoded with D = %d\n', D);
figure(2)
plot(x1(:,1),x1(:,2),'g.')
hold on
plot(x2(:,1),x2(:,2),'r*')
legend('Original','Codebook')
f = gcf;
exportgraphics(f,[savepic '2.png'])


%% Decode
x_end = int16(zeros(Nx,1)); 

idxT = 1;
for i = 1:L:Nx-L
    for j = 1:L
         x_end(i+j-1) = x2(idxT,j);
    end
    idxT = idxT + 1;
end

%% Compare 

figure(3)
plot(abs(x_end(1:Nx)-x))
f = gcf;
exportgraphics(f,[savepic '3.png'])

% sound(single(x)/(2^nbit),F); % original signal
figure(4);
plot(x(150000:160000,:),'g-');
hold on
% disp('Press any key to continue');
% pause;
% clear sound % Added to stop sound playing

% sound(single(x_end)/(2^nbit),F); % encoded signal
plot(x_end(150000:160000,:),'r-');
% disp('Press any key to continue');
% pause;
% clear sound % Added to stop sound playing
legend('Original','Encoded')
f = gcf;
exportgraphics(f,[savepic '4.png'])
disp('Saving Plots DONE')

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
                dist = norm(T(i,:)-y(j,:))^2;
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
            else                
                y(i,:) = sum(T(address==i,:))./cards(i);
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

function [x2,address, D] = quantizer(x1,b,L,K)
    N = length(x1(:,1));
    x2 = zeros(N,L);
    address = int16(zeros(N,1));
    cards = int16(zeros(K,1));
    D = 0;
    for i=1:N
        if mod(i,round(N/4)) == 0 % only show 4 times
            fprintf('Quantizier on index %d/%d\n',i,N)
        end
        mindist = inf;
        min_idx = 1;
        for j = 1:K
            dist = norm((x1(i,:)-b(j,:)))^2;
            if dist < mindist
                mindist = dist;
                min_idx = j;
            end
        end
        cards(min_idx) = cards(min_idx) + 1;
        address(i) = min_idx;
        x2(i,:) = b(min_idx,:);
        D = D + mindist;
    end
    D = D/N/L;      
end

function [x,F,Nx,maxX] = loadaudio(Naud,file)
    if nargin < 2 
        filename = ['Audio\' num2str(Naud)  'mono.wav'];
    else
        filename = file;
    end
    info = audioinfo(filename);
    [x,F] = audioread(filename,'native') ; 
%     fprintf('Sampling frequency:  F = %d [Hz] \n',F); 
%     fprintf('Resolution:          nbits = %d [bit] \n',info.BitsPerSample);
    
    % Upscale to 16 bit/sample
    if info.BitsPerSample == 8
%         disp('Upscaled to 16 bit/sample')
        x = int16(x) - 127;
    end
    
    % Convert to mono
    if info.NumChannels == 2
%         disp('Converted to Mono')
        x = int16(mean(x,2));
    end
    Nx = length(x);
    maxX = max(x);
end

function [x,F,Nx,maxX]  = loadAllAudio()
    Nauds = [54, 58, 64, 70];
    [x,F,~,~] = loadaudio(49);
    for Naud=Nauds
        [xi,~,~,~] = loadaudio(Naud);
        x = [x; xi];
    end
    Nx = length(x);
    maxX = max(x);
end

function [x,F,Nx,maxX]  = loadAllAudioMusic()
    [x,F,~,~] = loadAllAudio();
    musics = {'music\SayNada.wav'; 'music\Good4U.wav'; 'music\WaitingOnAWar.wav'; 
     'music\ForeverAfterAll.wav'; 'music\SummerThing.wav'; 'music\TodoDeTi.wav'; 'music\TooOfficial.wav'};
    for i=1:7
        [xi,~,~,~] = loadaudio(1,musics{i});
        x = [x; xi];
    end
    Nx = length(x);
    maxX = max(x);
end
