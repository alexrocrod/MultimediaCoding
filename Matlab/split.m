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