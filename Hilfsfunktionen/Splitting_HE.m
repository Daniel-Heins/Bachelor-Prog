%Split Lösung zu der Wärme Leitungs Gleichung mit Zeitschritt tau
function Split_first = Splitting_HE(v_0,tau,N)

    coef  = func_to_fourrier(v_0,N);
    z     = [exp(-((pi*(0:N-1)).^2)*tau).*coef zeros(1,N)];
    Split_first    = fourrier_to_func(z,N);

end
