function [x,F,Nx,maxX]  = loadAllAudio(L,K)
    Nauds = [54, 58, 64, 70];
    [x,F,~,~] = loadaudio(49,'',L,K/4);
    for Naud=Nauds
        [xi,~,~,~] = loadaudio(Naud,'',L,K/4);
        x = [x; xi];
    end
    Nx = length(x);
    maxX = max(x);
end