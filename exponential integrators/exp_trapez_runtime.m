%% code for figure 4.4

%% declaring functions for the example

func_u    = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_t  = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_xx = @(t, x)  60 * exp(t) .* x .* (1-x) .* (5*x.^2 - 5*x +1);  

%% g(t,x) = u_t(t,x)-u_xx(t,x), so u solves the original Problem

func_g    = @(t, x)   func_u_t(t,x) - func_u_xx(t,x);  


%% sandbox for variables       

time_step = 0.001;
T = 1;
L_vals = 7:12;

%% pre-allocating

    runtime_FFT = zeros(1,length(L_vals));
    runtime_BF  = zeros(1,length(L_vals));

%% time complexity of exponential trapez method as a function of L with FFT

p = 1;
for L = L_vals


    % setting of variables   
    N  = 2^L;
    x  = 0:1/N:(N-1)/N;              
    u_approx = func_u(0, x);
    exponent = exp(- time_step * (pi * (0:N-1)).^2);
    fourcoef_u = imag(fft (func_u(0, x) ));
    substit   = time_step * (pi*(1:N-1)).^2;

    % sequential splitting until T  
    tic;
    for t_old = 0:round(T/time_step)
        t_new = t_old + time_step;

        % fourriercoefficient of g until at t_old and t_new
            fourcoef_g_t_old = imag(fft(func_g(t_old,x)));
            fourcoef_g_t_new = imag(fft(func_g(t_new,x)));

        % approximation of the first term
            first_term  = exp(-((pi*(0:N-1)).^2)*time_step).* fourcoef_u;

        % approximation of the second term 
            b_2 = (exp(-substit) - 1 + substit) ./ (substit).^2;
            b_1 = ((1 - exp(-substit).^2)) ./ substit - b_2;
            second_term = fourcoef_g_t_old(1,2:end) .* b_1 + fourcoef_g_t_new(1,2:end) .* b_2;

        % approximation at t_new
            fourcoef_u = first_term + time_step * [0 second_term];
            u_approx  = imag(ifft(fourcoef_u,N));

    end
    runtime_FFT(p) = toc;
    p=p+1;
end
%% time complexity of exponential trapez method as a function of L with Brute-Force

p = 1;
for L = L_vals

    % setting of variables   
    N  = 2^L;
    x  = 0:1/N:(N-1)/N;              
    u_approx = func_u(0,x);
    exponent = exp( - time_step * ( pi * (0:N-1) ).^2 ) ;
    fourcoef_u = imag(fft (func_u(0, x) ));

    % Brute-Force Matrix
    n = 0:N-1;
    SinMatrix = pi * transpose(n) * n / N; 
    substit   = time_step * (pi*(1:N-1)).^2;
    % sequential splitting until T
    tic
    for t_old = 0:round(T/time_step)
        t_new = t_old + time_step;

        % fourriercoefficient of g until at t_old and t_new
            fourcoef_g_t_old = ( sqrt(2) / N ) * func_g(t_old,x) * SinMatrix;
            fourcoef_g_t_new = ( sqrt(2) / N ) * func_g(t_new,x) * SinMatrix;

        % approximation of the first term
            first_term  = exp(-((pi*(0:N-1)).^2)*time_step).* fourcoef_u;

        % approximation of the second term  
            b_2 = (exp(-substit) - 1 + substit) ./ (substit).^2;
            b_1 = ((1 - exp(-substit).^2)) ./ substit - b_2;
            second_term = fourcoef_g_t_old(1,2:end) .* b_1 + fourcoef_g_t_new(1,2:end) .* b_2;

        % approximation at t_new
            fourcoef_u = first_term + time_step * [0 second_term];
            u_approx  = fourcoef_u * SinMatrix;

    end
    runtime_BF(p) = toc;
    p=p+1;
end

%% Plotting

figure("Name","Runtime of exponential euler using FFT and brute-force", ...
    "NumberTitle","off","Position",[100, 100, 800, 600]);
plot(L_vals,runtime_FFT,"b-o","LineWidth",1.5);
hold on;
plot(L_vals,runtime_BF,"r-o","LineWidth",1.5);
hold off;
grid on;
set(gca ,"YScale" ,"log");
ylim([1e-2 1e4]);
xlabel("Potenz L");
ylabel("Laufzeiten (in Sek., logarithmisch)");
legend("Exp Trapez-Quad. mit FFT","Exp. Trapez-Quad. mit B.F.","Location","northwest");