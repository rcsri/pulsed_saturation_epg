function str = get_str_x(x_in, x_label)

str = [];
for ind = 1:length(x_label)
    
    if iscell(x_label)
        lab = x_label{ind};
    else
        lab = x_label(ind);
    end
    if strcmp(lab,'T1a')
        curr_str = sprintf('T1a(ms)=%5.1f', 1000/x_in(ind));
        
    elseif strcmp(lab,'T2a')
        curr_str = sprintf('T2a(ms)=%5.1f', 1000/x_in(ind));
        
    elseif strcmp(lab,'T1s')
        curr_str = sprintf('T1s(ms)=%5.1f', 1000/x_in(ind));
        
    elseif strcmp(lab,'T2s')
        curr_str = sprintf('T2s(us)=%3.1f', 1e6/x_in(ind));
        
    elseif strcmp(lab,'M0s')
        curr_str = sprintf('M0s=%5.4f', x_in(ind));
        
    elseif strcmp(lab,'Rs')
        curr_str = sprintf('Rs=%5.1f', x_in(ind));
        
    elseif strcmp(lab,'M0a')
        curr_str = sprintf('M0a=%3.1f', x_in(ind));
        
    elseif strcmp(lab,'T1c')
        curr_str = sprintf('T1c(ms)=%5.1f', 1000/x_in(ind));
        
    elseif strcmp(lab,'T2c')
        curr_str = sprintf('T2c(ms)=%5.1f', 1000/x_in(ind));
        
    elseif strcmp(lab,'M0c')
        curr_str = sprintf('M0c=%5.4f', x_in(ind));
        
    elseif strcmp(lab,'Rc')
        curr_str = sprintf('Rc=%5.1f', x_in(ind));
    elseif strcmp(lab,'B0')
        curr_str = sprintf('B0=%3.1f', x_in(ind));
    elseif strcmp(lab,'B1')
        curr_str = sprintf('B1=%3.1f', x_in(ind));
    else
        curr_str = 'na';
    end
    
    if ind > 1
        str = [str ', '];
    end
    str = [str curr_str];
end

