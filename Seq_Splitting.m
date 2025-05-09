function u = Seq_Splitting(Randfunktion, func_g, Zeitintervall, Diskret_in_x, Diskret_in_Zeit)


% Vereinfachung

    func_u0 = Randfunktion;
    T   = Zeitintervall;
    L   = Diskret_in_x;
    tau = Diskret_in_Zeit;

    
%Parametersetzung

    N = 2^L;
    x = 0:1/N:(N-1)/N;
    u(1,:) = func_u0(x);
    time_steps = round(T/tau);

    
%Splitting-Verfahren

    for n = 1:time_steps
        u_split_HE  = Splitting_HE (u(n,:),tau,N);
        u(n+1,:)    = Splitting_RK2(u_split_HE,func_g,tau);
    end

%Plotten

    t = 0:tau:(n)*tau;
    [X,Y] = meshgrid(x,t);
    figure("Name","Splitting");
    surf(X,Y,u);
    shading interp

clear ans;
end
