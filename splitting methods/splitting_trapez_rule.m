%Trapez Quadratur regel
function Split_second = splitting_trapez_rule(func_vars, func_g, time_step, t_old, x)
    
    t_new = t_old + time_step;
    Split_second = func_vars + (time_step/2) * (func_g(t_old, x)+func_g(t_new, x));

end

