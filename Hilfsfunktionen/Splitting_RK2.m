%RK Verfahren der Ordnung 2


function Split_second = Splitting_RK2(v,g,tau)


    k1 = g(v);
    k2 = g(v + tau * k1);
    
    Split_second = v + (tau/2) * (k1 + k2);
    
    
end
