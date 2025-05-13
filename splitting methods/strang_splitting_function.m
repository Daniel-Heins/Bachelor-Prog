function u = strang_splitting_function(Randfunktion, func_g, Zeitintervall, Diskret_in_x, Diskret_in_Zeit)
%% notation

    func_u0 = Randfunktion;
    T   = Zeitintervall;
    L   = Diskret_in_x;
    tau = Diskret_in_Zeit;

%% setting variables

    N = 2^L;
    x = 0:1/N:(N-1)/N;
    u(1,:) = func_u0(x);
    time_steps = round(T/tau);

%% marchuk-strang splitting method

    for n = 1:time_steps
        u_split_HE = Splitting_HE (u(n,:),tau/2,N);
        u_half     = Splitting_RK2(u_split_HE,func_g,tau);
        u(n+1,:)   = Splitting_HE (u_half,tau/2,N);
    end

%% Plotting

    t = 0:tau:n*tau;
    [X,Y] = meshgrid(x,t);
    figure("Name","Marchuk-Strang Splitting","NumberTitle","off");
    surf(X,Y,u);
    shading interp
    
    
end
