%% code for figure 4.3

%% declaring functions for the example

func_u    = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_t  = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_xx = @(t, x)  60 * exp(t) .* x .* (1-x) .* (5*x.^2 - 5*x +1);  

%% g(t,x) = u_t(t,x)-u_xx(t,x), so u solves the original Problem

func_g    = @(t, x)   func_u_t(t,x) - func_u_xx(t,x);  


%% sandbox for variables       

time_step = 0.001;
T = 1;
L_vals = 7:10;

%% pre-allocating

    runtime_FFT = zeros(1,length(L_vals));
    runtime_BF  = zeros(1,length(L_vals));

%% time complexity of exponential euler method as a function of L with FFT

p = 1;
for L = L_vals


    % setting of variables   
    N  = 2^L;
    x  = 0:1/N:(N-1)/N;              
    u_approx = func_u(0, x);
    exponent = exp(- time_step * (pi * (0:N-1)).^2);
    fourcoef_u = imag(fft (func_u(0, x) ));

    % sequential splitting until T  
    tic;
    for t_old = 0:round(T/time_step)
        t_new = t_old + time_step;

        % fourriercoefficient of g until at t_old
            fourcoef_g  = imag(fft (func_g(t_old, x)));

        % approximation of the first term
            first_term  = exponent .* fourcoef_u;

        % approximation of the second term 
            second_term = fourcoef_g(1, 2:end) .* exponent(1, 2:end) ./ (pi*(1:N-1)).^2;

        % approximation at t_new
            fourcoef_u  = first_term + [0 second_term];
            u_approx    = imag(ifft (fourcoef_u, N));
  

    end
    runtime_FFT(p) = toc;
    p=p+1;
end
%% time complexity of exponential euler method as a function of L with Brute-Force

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
    
    % sequential splitting until T
    tic
    for t_old = 0:round(T/time_step)
        t_new = t_old + time_step;

        % fourriercoefficient of g until at t_old
            fourcoef_g  = ( sqrt(2) / N ) * func_g(t_old, x) * SinMatrix;

        % approximation of the first term
            first_term  = exponent .* fourcoef_u;

        % approximation of the second term
            second_term = fourcoef_g(1, 2:end) .* (1 - exponent(1, 2:end)) ./ (pi * (1:N-1)).^2;

        % approximation at t_new
            fourcoef_u  = first_term + [0 second_term];
            u_approx    = fourcoef_u * SinMatrix;

    end
    runtime_BF(p) = toc;
    %Globaler Fehler über alle x_j, t_k.
    p=p+1;
end

%% Plotting

figure("Name","Runtime of exponential trapez using FFT and brute-force", ...
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
legend("Exp. Euler mit FFT","Strang Splitting mit B.-F.","Location","northwest");