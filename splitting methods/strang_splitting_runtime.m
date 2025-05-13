%% code for figure 3.4

%% sandbox for variables       

time_step = 0.001;
T = 1;
L_vals = 7:12;
%% declaring functions for the example

func_u    = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_t  = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_xx = @(t, x)  60 * exp(t) .* x .* (1-x) .* (5*x.^2 - 5*x +1);  

%% g(t,x) = u_t(t,x)-u_xx(t,x), so u solves the original Problem

func_g    = @(t, x)   func_u_t(t,x) - func_u_xx(t,x);  

%% pre-allocating

    runtime_FFT = zeros(1,length(L_vals));
    runtime_BF  = zeros(1,length(L_vals));

%% time complexity of marchuk-strang splitting as a function of L with FFT

p = 1;
for L = L_vals


    % setting of variables   
    N  = 2^L;
    x  = 0:1/N:(N-1)/N;              
    u_approx = func_u(0, x);
    exponent = - (pi * (0:N-1)).^2;

    % sequential splitting until T  
    tic;
    for t_old = 0:round(T/time_step)
        t_new = t_old + time_step;

        % seq. splitting for one time step
        fourrier_u   = imag(fft (u_approx, N));
        problem_A    = imag(fft (exp( t_new * exponent ) .* fourrier_u, N));

        problem_B    = Splitting_Trapez (problem_A, func_g, time_step, t_old, x);

        fourrier_u   = imag(fft (u_approx, N));
        problem_A    = imag(fft (exp( t_new * exponent ) .* fourrier_u, N));

        u_approx     = problem_B;
       
    end
    % save the running time
    runtime_FFT(p) = toc;
    p = p+1 ;
end

%% time complexity of marchuk-strang splitting as a function of L with Brute-Force

p = 1;
for L = L_vals

    % setting of variables   
    N  = 2^L;
    x  = 0:1/N:(N-1)/N;              
    u_approx = func_u(0,x);
    exponent = - (pi * (0:N-1)).^2 ;

    % Brute-Force Matrix
    n = 0:N-1;
    SinMatrix = pi * transpose(n) * n / N; 
    
    % sequential splitting until T
    tic
    for t_old = 0:round(T/time_step)
        t_new = t_old + time_step;

        % sequential splitting for one time step
        fourrier_u = ( sqrt(2) / N ) * u_approx * SinMatrix;
        problem_A  = ( exp( t_new * exponent ) .* fourrier_u ) * SinMatrix;

        problem_B  = Splitting_Trapez (problem_A, func_g, time_step, t_old, x);

        fourrier_u = ( sqrt(2) / N ) * problem_B * SinMatrix;
        problem_A  = ( exp( t_new * exponent ) .* fourrier_u ) * SinMatrix;

        u_approx  = problem_B;
       
    end
    % safe the running time
    runtime_BF(p) = toc;
    p = p+1 ;
end

%% Plotting 

figure("Name","Runtime of marchuk-strang splitting using FFT and brute-force", ...
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
legend("Strang Splitting FFT","Strang Splitting B.F.","Location","northwest");