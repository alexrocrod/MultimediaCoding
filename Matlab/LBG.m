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