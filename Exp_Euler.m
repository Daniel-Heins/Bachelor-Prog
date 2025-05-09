function u_app = Exp_Euler(Randfunktion, func_g, Zeitintervall, Diskret_in_x, Diskret_in_Zeit)

% Vereinfachung

    func_u0 = Randfunktion;
    T   = Zeitintervall;
    L   = Diskret_in_x;
    tau = Diskret_in_Zeit;

%Parametersetzung

    N = 2^L;
    x = 0:1/N:(N-1)/N;
    u_app(1,:) = func_u0(x);
    fourcoef_u = func_to_fourrier(u_app(1,:),N);
    t_new = 0;
    time_steps = round(T/tau);

%Exponentielle Euler Verfahren

    for n = 2 : time_steps
        %Zeitpunkte setzen
        t_old  = t_new;
        t_new  = t_old + tau;

        %Fourrierkoeffizienten von g zur Zeit t_old bestimmen
            fourcoef_gu = func_to_fourrier(func_g(u_app(n-1,:)),N);

        %Approximation von exp(-hA)u^{k+1}
            first  = exp(-((pi*(0:N-1)).^2)*tau).* fourcoef_u;
            
        %Approximation des Integrals der Variation der Konstanten Formel 
            second = fourcoef_gu(1,2:end) .* (1 - exp(-tau*(pi*(1:N-1)).^2)) ./ (pi*(1:N-1)).^2;
            
        %Approximation der Lösung zum Zeitpunkt t_new
            fourcoef_u = first + [0 second];
            
        %Exakte Werte für Vergleich berechnen
            u_app(n,:)  = fourrier_to_func(fourcoef_u,N);
            
    end

%Plotten

    t = 0:tau:(n-1)*tau;
    [X,Y] = meshgrid(x,t);
    figure("Name","Exponentielle Euler Verfahren","NumberTitle","off");
    surf(X,Y,u_app);
    shading interp

end
