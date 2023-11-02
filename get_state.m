function state = get_state(Pop, dlow, dhigh)
    % 获取当前种群的状态
    Xmax = 100;
    Xmin = -100;
%     
    Pop = (Pop-Xmin)/(Xmax - Xmin);
    center_pop = mean(Pop,1);
    
    ps = size(Pop, 1);
    div = 0.0;

    
    for i=1:ps
        div = div + sqrt(sum((Pop(i,:)-center_pop).^2));

    end

    diversity = div/(ps);
    if diversity > dhigh
        state = 1;
    elseif diversity < dlow
        state = 3;
    else
        state = 2;
    end
end