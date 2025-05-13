%% code for figure 3.1
clear
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

%% time complexity of sequential splitting as a function of L with FFT
p = 1;
for time_step = Steps

    % pre-allocating
    error_at_t_new = zeros(1, round(T/time_step));

    % sequential splitting until T 
    k = 1;
    t_new = 0;
    u_approx = func_u(0, x);

    for index = 1:round(T/time_step)
        t_old = t_new;
        t_new = t_old + time_step;

        % seq. splitting for one time step
        problem_A    = Splitting_HE (u_approx, time_step, N_val);
        problem_B    = Splitting_Trapez (problem_A, func_g, time_step, t_old, x);

        u_approx     = problem_B;

        % safe the error at t_new
        error_at_t_new(index)   = norm(u_approx - func_u(t_new, x), inf) ;
    end
    %globaler Fehler 
    global_error(p) = norm(error_at_t_new,inf);
    p = p+1 ;
end

%% Plotten
 
figure("Name","Seq. Splitting Fehler", ...
    "NumberTitle","off","Position",[100, 100, 800, 600]);
% constants for nice plots
    C_1 = 10;      
    C_2 = 0.01;  
obere_schranke     = C_1 * (Steps) ; 
untere_schranke    = C_2 * (Steps.^2);
loglog(Steps, obere_schranke, "g--", "Linewidth", 1.5);
hold on;
loglog(Steps, global_error, "r-o", "Linewidth", 2);
hold on;
loglog(Steps, untere_schranke, "b--", "Linewidth", 1.5);
grid on;
ylim([1e-8 1e1]);
xlabel("Schrittweiten (\tau, logarithmisch)");
ylabel("globaler Fehler (logarithmisch)");
legend("C_1 tau  ,C_1 ∊ ℝ " ...
    ,"Seq. Splitting Fehler","C_2 tau^2  ,C_2 ∊ ℝ ","Location","northwest");