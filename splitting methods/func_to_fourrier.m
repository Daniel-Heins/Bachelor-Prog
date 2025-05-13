%% turn function to fourriercoefficients

function coef = func_to_fourrier(funk_werte,N)

    z = [funk_werte zeros(1,N)];
    FFT_z = sqrt(2)*imag(FFt(z,2*N))/N;
    coef    = FFT_z(1,1:N);
    
end

