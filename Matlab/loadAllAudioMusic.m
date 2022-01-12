function [x,F,Nx,maxX]  = loadAllAudioMusic(L,K)
    [x,F,~,~] = loadAllAudio(L,round(K*4/11));
    musics = {'music\SayNada.wav'; 'music\Good4U.wav'; 'music\WaitingOnAWar.wav'; 
     'music\ForeverAfterAll.wav'; 'music\SummerThing.wav'; 'music\TodoDeTi.wav'; 'music\TooOfficial.wav'};
    for i=1:7
        [xi,~,~,~] = loadaudio(1,musics{i},L,round(K*7/11));
        x = [x; xi];
    end
    Nx = length(x);
    maxX = max(x);
end