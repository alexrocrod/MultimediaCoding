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
playsound = false; % if true both final signals are played
epsilon = 1e-2;
L = 2; % 2 or 4
iL = 1; % index to identify the rate selected
Naud = 100; % Training Set: 100- All Audio, 200- SayNada, 300- All Audio and Music
NEnc = 49; % Encoding Object ID: [49,54,58,64,70] Audio, 200- SayNada,
           % 210- Good4U, 220- WaitingOnAWar, 230- ForeverAfterAll, 
           % 240- SummerThing, 250- TodoDeTi, 260- TooOfficial

if L == 2, Rs = [2, 4];
elseif L == 4, Rs = [1, 2];
else, disp('Invalid vector dimension L, valid values are 2 or 4.'); return;
end

Ks = 2.^(Rs.*L); % Codebook sizes

nbit = 16;
K = Ks(iL);
R = Rs(iL);

% Saving Codebooks and Results location
savefile = ['CodeBooks\cb_' num2str(L)  '_' num2str(R) '_' num2str(epsilon) '_' num2str(Naud) '.mat'];
savepic = ['Results\cb_' num2str(L)  '_' num2str(R) '_' num2str(epsilon) '_' num2str(Naud) '\Enc' num2str(NEnc) '\'];
mkdir(savepic(1:end-1))

%% Load Training Audio

if Naud == 100 % All Audio
    [x,F,Nx,maxX] = loadAllAudio(L,K);
elseif Naud == 200 % Say Nada
    [x,F,Nx,maxX] = loadaudio(1,'music\SayNada.wav',L,K);
elseif Naud == 300 % All audio and Music
    [x,F,Nx,maxX] = loadAllAudioMusic(L,K);
else
    disp('Invalid Training Set ID, valid values: 100, 200, 300.')
    return
end

%% Train

if isfile(savefile)
    load(savefile)
    fprintf('Codebook loaded from file %s\n',savefile)
else
    % Training set 
    N = round(Nx/L)-1;
    T = zeros(N,L); 
    
    idxT = 1;
    for i = 1:L:Nx-L
        for j = 1:L
            T(idxT,j) = x(i+j-1); 
        end
        idxT = idxT + 1;
    end
    fprintf('N/K = %d\n', N/K);

    % Compute Codebook
    tic
    y = split(T,epsilon,K,maxX,L);
    toc 
    save(savefile,'y')
    
    % Plot Codebook over Original Vectorized Signal
    plot(T(:,1), T(:,2), 'g.', y(:,1), y(:,2), 'r*')
    legend('Original','Codebook')
    f = gcf;
    exportgraphics(f, [savepic '1.png'])
end

%% Load Encoding Audio
if NEnc <= 70, [x,F,Nx,maxX] = loadaudio(NEnc); % Audio
elseif NEnc == 200, [x,F,Nx,~] = loadaudio(1,'music\SayNada.wav');
elseif NEnc == 210, [x,F,Nx,~] = loadaudio(1,'music\Good4U.wav');
elseif NEnc == 220, [x,F,Nx,~] = loadaudio(1,'music\WaitingOnAWar.wav');
elseif NEnc == 230, [x,F,Nx,~] = loadaudio(1,'music\ForeverAfterAll.wav');
elseif NEnc == 240, [x,F,Nx,~] = loadaudio(1,'music\SummerThing.wav');
elseif NEnc == 250, [x,F,Nx,~] = loadaudio(1,'music\TodoDeTi.wav');
elseif NEnc == 260, [x,F,Nx,~] = loadaudio(1,'music\TooOfficial.wav');
end
fprintf('Variance of x = %d\n', var(double(x)));


%% Encode

% Read vectors from original signal
tic
x1 = zeros(round(Nx/L),L); % vectorized x before quantization
idxT = 1;
for i = 1:L:Nx-L
    for j = 1:L
        x1(idxT,j) = x(i+j-1); 
    end
    idxT = idxT + 1;
end

% Quantize audio signal
[x2,ads,D] = quantizer(x1,y,L,K); % x2- vectorized x after quantization
toc % Encoding Time

% Plot Quantized over Original 
fprintf('Audio Encoded with D = %d\n', D);
figure(2)
plot(x1(:,1),x1(:,2),'g.', x2(:,1),x2(:,2),'r*')
legend('Original','Encoded')
f = gcf;
exportgraphics(f,[savepic '2.png'])


%% Decode
% From vector form to Mono Audio Signal
x_end = int16(zeros(Nx,1)); 
idxT = 1;
for i = 1:L:Nx-L
    for j = 1:L
         x_end(i+j-1) = x2(idxT,j);
    end
    idxT = idxT + 1;
end

%% Compare 
% Plot absolute error for each data point
figure(3)
plot(abs(x_end(1:Nx)-x))
f = gcf;
exportgraphics(f,[savepic '3.png'])

% Plot Original Signal
figure(4);
plot(x(150000:160000,:),'g-');
hold on
if playsound % Play Original
    sound(single(x)/(2^nbit),F); % original signal
    disp('Press any key to continue');
    pause;
    clear sound % Added to stop sound playing
end

% Plot Final Signal
plot(x_end(150000:160000,:),'r-');
if playsound % Play Final
    sound(single(x_end)/(2^nbit),F); % encoded signal
    disp('Press any key to continue');
    pause;
    clear sound % Added to stop sound playing
end

legend('Original','Encoded')
f = gcf;
exportgraphics(f,[savepic '4.png'])
