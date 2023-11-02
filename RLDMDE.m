function [gbestval,gbest,BestCost] = RLDMDE( fhd,Particle_Number,Dimension,Popmin,Popmax,Max_nfe,varargin)
    warning('off');
%     rand('seed', sum(100 * clock));

    %% 变量定义
    ps = Particle_Number;  
    D = Dimension; 
    Gen = round(Max_nfe/ps); 
    K = 3; 
    ps_K = ps/K; %
    VarMin= Popmin;          % Lower Bound of Decision Variables  
    VarMax= Popmax;          % Upper Bound of Decision Variables
    FES = 0; 
    
    d_low = 5e-6; 
    d_high = 0.25; 
    
    state_num = 3;
    action_num = 3;
    gamma = 0.8; 
    alpha = 0.2; 
    
%     BestCost = zeros(1, Gen);
    
   %% 初始化
    m_pos = Pos_init(ps,D,VarMax,VarMin);
    m_opos = VarMin+VarMax-rand(ps,D).*m_pos;
    m_pos = [m_pos;m_opos];
    m_eval = feval(fhd,m_pos',varargin{:}); % 种群评估
    
    [~,index] = sort(m_eval);
    m_pos = m_pos(index(1:ps),:);
    m_eval = m_eval(index(1:ps));
    m_eval = m_eval';
    
    
    del = randperm(ps);
    m_pos = m_pos(del,:);
    m_eval = m_eval(del);
    FES = 2*ps; % 更新迭代次数
    [gbestval,gbestid] = min(m_eval);
    gbest = m_pos(gbestid,:);
    BestCost(1) = gbestval; 
    
    

    %% 小组定义
    state_max = K;
    
    Pop=repmat(struct('Pos',{},'eval',{},'best',{},'b_val',{},'Qtable',{},'state',{},'state_record',{},'cnt_r',{},'ps_num',{}),1,K);
    for i=1:K
        Pop(i).Pos = m_pos(ps_K*(i-1)+1:ps_K*i, :);
        Pop(i).eval = m_eval(ps_K*(i-1)+1:ps_K*i);
        [Pop(i).b_val,id] = min(Pop(i).eval);
        Pop(i).best = Pop(i).Pos(id,:);
       
        
        Pop(i).state =  get_state(Pop(i).Pos, d_low, d_high);
        Pop(i).Qtable = zeros(state_num,action_num);
        Pop(i).state_record = zeros(1,state_max);
        Pop(i).cnt_r = 1;
        Pop(i).state_record(Pop(i).cnt_r) = Pop(i).state;
        Pop(i).ps_num = ps_K;
    end
    
   %% SHADE　参数设置
    memory_size = ps;
    memory_sf = 0.5 .* ones(memory_size, 1);
    memory_cr = 0.5 .* ones(memory_size, 1);
    memory_pos = 1;
    
    
    %% 外部存储
    arc_rate = 2;
    archive.NP = arc_rate*ps; % the maximum size of the archive
    archive.pop = zeros(0, D); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions
    
    %% 主循环
    iter = 1;
    
    state_flag = zeros(1,K);
    stag_flag = zeros(1,K);
    
    
   
    
    while FES <= Max_nfe
        iter = iter + 1;
        mem_rand_index = ceil(memory_size * rand(ps, 1));
        mu_sf = memory_sf(mem_rand_index);
        mu_cr = memory_cr(mem_rand_index);
        
        %% for generating crossover rate
          cr = normrnd(mu_cr, 0.1);
          term_pos = find(mu_cr == -1);
          cr(term_pos) = 0;
          cr = min(cr, 1);
          cr = max(cr, 0);

          %% for generating scaling factor
          sf = mu_sf + 0.1 * tan(pi * (rand(ps, 1) - 0.5));
          pos = find(sf <= 0);

          while ~ isempty(pos)
            sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
            pos = find(sf <= 0);
          end

          sf = min(sf, 1); 
          
          
          %% 定义成功参加变异个体的F和CR集合
            goodF = [];
            goodCR = [];
            dif_val = [];
            
            uPos = [];
            ueval = [];
          
          %% 小组更新
          cnt = 0;
           
          ps_record = zeros(1,K);
          for i=1:K
              
              ps_record(i) = Pop(i).ps_num;
              if Pop(i).ps_num == 0
                  continue;
              end
              F = sf(cnt+1:cnt+Pop(i).ps_num);
              CR = cr(cnt+1:cnt+Pop(i).ps_num);
              popAll = [ Pop(i).Pos; archive.pop];
              action = choose_Action(Pop(i).Qtable, Pop(i).state);
              Q_predict = Pop(i).Qtable(Pop(i).state, action); 
              
%               disp(size(Pop(i).Pos));
              newPop = perform_action(action, Pop(i).Pos, Pop(i).eval, popAll, F, CR, VarMax, VarMin);
              neweval = feval(fhd,newPop',varargin{:});
              neweval = neweval';
              
              for j=1:Pop(i).ps_num
                  FES = FES + 1;
                  if neweval(j) < gbestval
                      gbest = newPop(j,:);
                      gbestval = neweval(j);
                  end
                  if FES > Max_nfe
                      break;
                  end
              end
              
              %% I == 1: the parent is better; I == 2: the offspring is better
                dif = abs(Pop(i).eval - neweval);
                I = (Pop(i).eval > neweval);
                
                
                goodCR = [goodCR;cr(I == 1)];  
                goodF = [goodF;sf(I == 1)];
                dif_val = [dif_val;dif(I == 1)];
               
                suc_rate = sum(I)/Pop(i).ps_num;
                
                uPos = [uPos; Pop(i).Pos(I == 1, :)];
                ueval = [ueval;Pop(i).eval(I == 1)];
                
                
                [Pop(i).eval, I] = min([Pop(i).eval, neweval], [], 2);
                Pop(i).Pos(I==2,:) = newPop(I==2,:);
                
                [nv,id] = min(Pop(i).eval);
                if nv < Pop(i).b_val
                    Pop(i).b_val = nv;
                    Pop(i).best = Pop(i).Pos(id,:);
                end
                 
                 if suc_rate <= 1/Pop(i).ps_num
                    R1 = -1*suc_rate;
                elseif suc_rate > 1/Pop(i).ps_num && suc_rate <= 0.5
                    R1 = 6*suc_rate;
                elseif suc_rate > 0.5 && suc_rate <= 0.75
                    R1 = 8*suc_rate;
                else
                    R1 = 10*suc_rate; % 将进化成功作为奖励
                 end
                
                S_n = get_state(Pop(i).Pos, d_low, d_high);
                Q_target = R1 + gamma * max(Pop(i).Qtable(Pop(i).state, :));
                Pop(i).Qtable(Pop(i).state, action) = Pop(i).Qtable(Pop(i).state, action) + alpha * (Q_target - Q_predict);
                Pop(i).state = S_n;
                cnt = cnt + Pop(i).ps_num;
                
                %% 状态记录
                if Pop(i).state ~= Pop(i).state_record(Pop(i).cnt_r) % 状态存在变化
                    Pop(i).cnt_r = 1;
                    Pop(i).state_record = zeros(1,state_max); %清空
                    Pop(i).state_record(Pop(i).cnt_r) = Pop(i).state;
                    state_flag(i) = 0;
                else
                    Pop(i).cnt_r = Pop(i).cnt_r + 1;
                    if Pop(i).cnt_r == state_max % 5次状态均没有变化
                        state_flag(i) = 1; % 存在停滞
                        Pop(i).state_record = zeros(1,state_max); %清空
                        Pop(i).cnt_r = 1;
                    else
                        Pop(i).state_record(Pop(i).cnt_r) = Pop(i).state;
                    end
                    
                end
                
                %% 停滞判断
                if state_flag(i) == 1 && suc_rate <= 1/Pop(i).ps_num
                    stag_flag(i) = 1;
                else
                    stag_flag(i) = 0;
                end
                
                
          end
%         disp(ps_record)
%         disp(stag_flag)
       %%  更新全局参数
       archive = updateArchive(archive, uPos, ueval);
       num_success_params = numel(goodCR);     
       if num_success_params > 0 
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;

        %% for updating the memory of scaling factor 
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);

        %% for updating the memory of crossover rate
        if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
          memory_cr(memory_pos)  = -1;
        else
          memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
        end

        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1; end
       end
      
       %% 停滞处理
       % 存在优势小组
       if sum(stag_flag) ~= K
          mig_counter = 0;
          add_pos = [];
          add_eval = [];
          %% 劣势小组个体迁出
          for i=1:K
              if stag_flag(i) == 1
                  newpopsize = floor(Pop(i).ps_num*0.9); %新的种群大小
                  ct = max(0,Pop(i).ps_num-newpopsize);
                  if ct == 0
                      continue
                  end
                  if newpopsize == 1
                      ct = 2; 
                      newpopsize = 0;
                      ps_record(i) = 0;
                  end
                  mig_counter = mig_counter + ct;
                  del_r = randperm(Pop(i).ps_num, ct);
                  add_pos = [add_pos; Pop(i).Pos(del_r,:)];
                  add_eval = [add_eval; Pop(i).eval(del_r)];
                  Pop(i).Pos(del_r,:) = [];
                  Pop(i).eval(del_r) = [];
                  Pop(i).ps_num = newpopsize;
                  
                  %% 产生粒子删减后，状态重新记录
                   Pop(i).state = get_state(Pop(i).Pos, d_low, d_high);
                   Pop(i).state_record(Pop(i).cnt_r) = Pop(i).state;
                   [Pop(i).b_val,id] = min(Pop(i).eval);
                   Pop(i).best = Pop(i).Pos(id,:);
%                    stag_flag(i) = 0; % 停滞清零
              end
          end
          
          %% 优势小组吸收劣势小组个体
          if mig_counter ~= 0 
              ps_All = sum(~stag_flag.*ps_record); %优势小组的全部个体
              
              c_flag = ~stag_flag;
              
              tt = mig_counter;
              for i=1:K

                  if stag_flag(i) == 1
                      continue
                  end
                  if sum(c_flag) == 1
                      % 表明只有一个组受到奖励
                      in_ps = tt;
                  else
                      in_ps = round(mig_counter*ps_record(i)/ps_All);
                  end
                  if in_ps <= 0
                      continue
                  end
                  
                  Pop(i).Pos = [Pop(i).Pos;add_pos(1:in_ps,:)];
                  Pop(i).eval = [Pop(i).eval;add_eval(1:in_ps)];
                  Pop(i).ps_num = size(Pop(i).Pos,1);
                  add_pos(1:in_ps,:) = [];
                  add_eval(1:in_ps) = [];
                  tt = tt - in_ps;
                  Pop(i).state = get_state(Pop(i).Pos, d_low, d_high);
                  Pop(i).state_record(Pop(i).cnt_r) = Pop(i).state;
                 [Pop(i).b_val,id] = min(Pop(i).eval);
                  Pop(i).best = Pop(i).Pos(id,:);
                  c_flag(i) = 0;
              end

          end
       
          
       else
%            %% 种群均存在停滞
           m_pos = [];
           m_eval = [];
           for i=1:K
               m_pos = [m_pos;Pop(i).Pos];
               m_eval = [m_eval;Pop(i).eval];
           end
           
           
%            
            del = randperm(ps);
            m_pos = m_pos(del,:);
            m_eval = m_eval(del);
            
            
%            c_flag = zeros(1,K); 
           
           for i=1:K
               
               if Pop(i).ps_num == 0
                   continue
               end
               
               
               
                Pop(i).Pos = m_pos(1:Pop(i).ps_num, :);
                Pop(i).eval = m_eval(1:Pop(i).ps_num);
                m_pos(1:Pop(i).ps_num, :) = [];
                m_eval(1:Pop(i).ps_num) = [];
                
                [Pop(i).b_val,id] = min(Pop(i).eval);
                Pop(i).best = Pop(i).Pos(id,:);
                Pop(i).state =  get_state(Pop(i).Pos, d_low, d_high);

                Pop(i).state_record = zeros(1,state_max);
                Pop(i).cnt_r = 1;
                Pop(i).state_record(Pop(i).cnt_r) = Pop(i).state;
                

            end
           
       end
       
        BestCost(iter) = gbestval;
    end

end