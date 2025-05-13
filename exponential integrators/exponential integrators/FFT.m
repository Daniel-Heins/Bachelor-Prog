%% FFT-Algorithm (Cooley and Tukey)

function Y=FFT(x,N)

w=exp(1i*2*pi/N);
Y=zeros(1,N);

if N==1
    Y=x;
else
    Y_u=FFt(x(2*(1:N/2)-1),N/2);
    Y_v=FFt(x(2*(1:N/2)),N/2);
    for k=1:N/2
        Y(k)=Y_u(k) + w^(k-1)*Y_v(k);
        Y(k+N/2)=Y_u(k) - w^(k-1)*Y_v(k);
    end
end

end
