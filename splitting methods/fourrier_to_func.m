%% helper function

%% turn fourriercoefficients to function

function func = fourrier_to_func(fourrier_werte,N)

    z = [fourrier_werte zeros(1,N)];
    FFT_z = sqrt(2)*imag(FFt(z,2*N));
    func = FFT_z(1,1:N);
    
end
