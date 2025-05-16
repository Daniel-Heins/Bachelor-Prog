%% figure 5.2

%% sandbox for variables
    T  = 10;
    x0 = 0.2;
    step_size  = 0.1;
    time_steps = round(T / step_size);

 %% pre-allocating
    x_vals = zeros(1,length(time_steps));
    x_exact = zeros(1,length(time_steps));

    t_new      = 0;

    x_vals(1)       = x0;
    x_exact(1)   = x0;
    
%% Solution of the ODE
    func_x = @(t) (x0 - (2/5)) * exp(-2*t) + (2/5) * cos(t) + (1/5) * sin(t); 

%% Sequential splitting method
    for n = 2 : time_steps+1
        t_old = t_new;
        t_new = t_new + step_size;
      
            x_A = exp(-2*step_size)*x_vals(n-1);
            x_B = x_A + sin(t_new) - sin(t_old);
    
            x_vals(n) = x_B;
            x_exact(n) = func_x(t_new);
    end

%% Plotting
    t_disk = 0:T/time_steps:T;
    figure("Name","Seq. Splitting Fehler","NumberTitle","off");
    plot(t_disk,x_vals,"r-","Linewidth",1.5);
    hold on;
    plot(t_disk,x_exact,"b-","Linewidth",1.5);
    grid on;
    ylim([-0.6 0.6]);
    xlabel("x-Werte");
    legend("exakte Lösung","App. Lösung mittels seq. Splitting","Location","northwest")
