%% code for figure 4.2

%% declaring functions for the example

func_u    = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_t  = @(t, x)  10 * exp(t) .* x.^3 .* (1-x).^3;
func_u_xx = @(t, x)  60 * exp(t) .* x .* (1-x) .* (5*x.^2 - 5*x +1);  

%% g(t,x) = u_t(t,x)-u_xx(t,x), so u solves the original Problem

func_g    = @(t, x)   func_u_t(t,x) - func_u_xx(t,x);  


%% sandbox for variables       

Steps = 2.^(0:7) * 0.001;
T = 1;
L_val = 9;     %% N = 256 

%% setting of variables   

N_val  = 2^L_val;
x  = 0 : 1/N_val : (N_val-1)/N_val;              
t_new = 0;
global_error = zeros(1, length(Steps));

%% time complexity of exponential trapez method as a function of L with FFT

p = 1;
for time_step = Steps

    % pre-allocating
    error_at_t_new = zeros(1, round(T/time_step)); 

    % setting variables
    t_new = 0;
    u_approx = func_u(0, x);
    substit   = time_step * (pi*(1:N_val-1)).^2;
    fourcoef_u = func_to_fourrier(func_u(0, x), N_val);

    %% trapez quadratur rule
    for index = 1:round(T/time_step)
        t_old = t_new;
        t_new = t_old + time_step;

        % fourriercoefficient of g until at t_old and t_new
            fourcoef_g_t_old = func_to_fourrier(func_g(t_old,x), N_val);
            fourcoef_g_t_new = func_to_fourrier(func_g(t_new,x), N_val);


        % approximation of the first term
            first_term  = exp(-((pi*(0:N_val-1)).^2)*time_step).* fourcoef_u;

        % approximation of the second term 
            b_2 = (exp(-substit) - 1 + substit) ./ (substit).^2;
            b_1 = ((1 - exp(-substit))) ./ substit - b_2;
            second_term = fourcoef_g_t_old(1,2:end) .* b_1 + fourcoef_g_t_new(1,2:end) .* b_2;

        % approximate the solution
            fourcoef_u = first_term + time_step * [0 second_term];

        % compare to solution
            u_approx  = fourrier_to_func(fourcoef_u, N_val);
            u_exakt   = func_u(t_new, x);

        % safe the error at t_new
            error_at_t_new(index) = norm(u_exakt-u_approx, inf);

    end
    % global error
        global_error(p) = norm(error_at_t_new,inf);
    p=p+1;
end

%% Plotten

figure("Name","Exp. Trapez Quadratur Fehler", ...
    "NumberTitle","off","Position",[100, 100, 800, 600]);
% constants for nice plot
    C_1 = 10;       
    C_2 = 0.01;    
obere_schranke     = C_1 * (Steps) ; 
untere_schranke    = C_2 * (Steps.^2); 
loglog(Steps,obere_schranke,"g--","Linewidth",1.5);
hold on;
loglog(Steps,global_error,"r-o","Linewidth",2);
hold on;
loglog(Steps,untere_schranke,"b--","Linewidth",1.5);
grid on;
ylim([1e-8 1e1]);
xlabel("Schrittweiten (\tau, logarithmisch)");
ylabel("globaler Fehler (logarithmisch)");
legend("C_1 tau^2  ,C_1 ∊ ℝ " ...
    ,"Exp. Trapez Fehler","C_2 tau  ,C_2 ∊ ℝ  ","Location","northwest");
