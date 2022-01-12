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