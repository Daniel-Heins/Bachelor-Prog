function u_app = exp_euler_function(Randfunktion, func_g, Zeitintervall, Diskret_in_x, Diskret_in_Zeit)
%% notations

    func_u0 = Randfunktion;
    T   = Zeitintervall;
    L   = Diskret_in_x;
    tau = Diskret_in_Zeit;

%% setting variables

    N = 2^L;
    x = 0:1/N:(N-1)/N;
    u_app(1,:) = func_u0(x);
    fourcoef_u = func_to_fourrier(u_app(1,:),N);
    t_new = 0;
    time_steps = round(T/tau);

%% exponential euler method

    for n = 2 : time_steps
        t_old  = t_new;
        t_new  = t_old + tau;

        % fourriercoefficients of g(u(t_old)) 
            fourcoef_g = func_to_fourrier(func_g(u_app(n-1,:)),N);

        % approximation of first term
            first_term  = exp(-((pi*(0:N-1)).^2)*tau).* fourcoef_u;
            
        % approximation of second term 
            second_term = fourcoef_g(1,2:end) .* (1 - exp(-tau*(pi*(1:N-1)).^2)) ./ (pi*(1:N-1)).^2;
            
        % approximation of fourriercoefficients at t_new
            fourcoef_u = first_term + [0 second_term];
            
        % approximation of solution at t_new
            u_app(n,:)  = fourrier_to_func(fourcoef_u,N);
            
    end

%% Plotting

    t = 0:tau:(n-1)*tau;
    [X,Y] = meshgrid(x,t);
    figure("Name","Exponentielle Euler Verfahren","NumberTitle","off");
    surf(X,Y,u_app);
    shading interp

end
