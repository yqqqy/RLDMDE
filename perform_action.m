function newPop = perform_action(action, Pop, val,  popAll, F, Cr, Xmax, Xmin)
    %% 该文件用于确定执行哪种突变策略
    p = 0.1;
    % 初始0.1
    [popsize, D] = size(Pop); %种群大小 维度
    rot = (0:1:popsize-1);
    
    
     if action == 1
         %%  DE/rand/1/bin
         ind = randperm(2);
         a1  = randperm(popsize);             % shuffle locations of vectors
         rt = rem(rot+ind(1),popsize);        % rotate indices by ind(1) positions
         a2  = a1(rt+1);                 % rotate vector locations
         rt = rem(rot+ind(2),popsize);
         a3  = a2(rt+1);
         posr1 = Pop(a1,:);             % shuffled population 1
         posr2 = Pop(a2,:);             % shuffled population 1
         posr3 = Pop(a3,:);             % shuffled population 1
         
%          gbestrep= repmat(gbest,popsize,1);
         vi = posr1 + F(:, ones(1,D)).*(posr2-posr3); % 突变向量
%          vi = gbestrep + F(:, ones(1,D)).*(posr2-posr3); % 突变向量
         
         vi = ((vi>=Xmin)&(vi<=Xmax)).*vi...
                +(vi<Xmin).*(Xmin+Pop)/2+...
                +(vi>Xmax).*(Xmax+Pop)/2; % 边界处理
         
        mask = rand(popsize, D) > Cr(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
        rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
        jrand = sub2ind([popsize D], rows, cols); mask(jrand) = false;
        newPop = vi; newPop(mask) = Pop(mask);
         
     elseif action == 2
         
         %% DE/current-to-pbest/1/bin
         r0 = [1 : popsize];
         [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
         [~, indBest] = sort(val, 'ascend');
         pNP = max(round(p * popsize), 2); % choose at least two best solutions
         randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
         randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
         pbest = Pop(indBest(randindex), :); % randomly choose one of the top 100p% solutions
         vi = Pop + F(:, ones(1, D)) .* (pbest - Pop + Pop(r1, :) - popAll(r2, :));
         vi = ((vi>=Xmin)&(vi<=Xmax)).*vi...
                +(vi<Xmin).*(Xmin+Pop)/2+...
                +(vi>Xmax).*(Xmax+Pop)/2; % 边界处理
        mask = rand(popsize, D) > Cr(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
        rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
        jrand = sub2ind([popsize D], rows, cols); mask(jrand) = false;
        newPop = vi; newPop(mask) = Pop(mask);
     else
         %%  DE/pbest/1/bin
         ind = randperm(1);
         a1  = randperm(popsize);             % shuffle locations of vectors
         rt = rem(rot+ind(1),popsize);        % rotate indices by ind(1) positions
         a2  = a1(rt+1);                 % rotate vector locations
         posr1 = Pop(a1,:);             % shuffled population 1
         posr2 = Pop(a2,:);             % shuffled population 2
        [~, indBest] = sort(val, 'ascend');
         pNP = max(round(p * popsize), 2); % choose at least two best solutions
         randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
         randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
         pbest = Pop(indBest(randindex), :); % randomly choose one of the top 100p% solutions
         
         vi = Pop + F(:, ones(1,D)).*(pbest-Pop) + rand(popsize,D).*(posr1-posr2); % 突变向量


         vi = ((vi>=Xmin)&(vi<=Xmax)).*vi...
                +(vi<Xmin).*(Xmin+Pop)/2+...
                +(vi>Xmax).*(Xmax+Pop)/2; % 边界处理
%          newPop = vi;
        mask = rand(popsize, D) > Cr(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
        rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
        jrand = sub2ind([popsize D], rows, cols); mask(jrand) = false;
        newPop = vi; newPop(mask) = Pop(mask);
     end

end