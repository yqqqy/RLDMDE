function action = choose_Action(Q_table, state)
    % ѡ�����Ŷ���
    epsilon = 0.8;

    r = rand;
    choices = Q_table(state, :);
    
    if r>epsilon || ~any(choices)
        action = randi(size(Q_table,2));
    else
        [~ ,action] = max(choices);
    end
    
   
end