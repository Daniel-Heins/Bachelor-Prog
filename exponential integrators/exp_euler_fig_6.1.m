%% figure 6.1

%% sandbox for variables
    A = [-2 0; 0 3];
    B = [0 2; -4 0];
    
    u0 = [1; 0];   
    T = 1;         
    
    step_size_array = 2.^(0:7) * 0.001;

%% pre-allocating
    errors = zeros(size(step_size_array)); 

%% Iteration for stept_size_array
    for k = 1:length(step_size_array)
        step_size = step_size_array(k);
        N = round(T / step_size);
        
    %% initialization
        u = zeros(2, N+1);
        u(:,1) = u0;
        t = linspace(0, T, N+1);
        
        expA = expm(A * step_size);
        phi1 = (expA - eye(size(A))) / (A * step_size);

    %% exponential euler method 
        for n = 1:N
            u(:,n+1) = expA * u(:,n) + step_size * phi1 * B * u(:,n);
        end
    
    %% exact solution
        u_exact = zeros(2, N+1);
        for n = 1:N+1
            u_exact(:,n) = expm((A+B)*t(n)) * u0;
        end
    
    %% global error 
        diff = u - u_exact;
        errors(k) = max(vecnorm(diff, inf)); 
    end

%% Plotting
    figure("Name","error plot exp. euler method for matrix", ...
        "NumberTitle","off"); 
    % constant for nice plot
    C_1 = 50;      
    untere_schranke    = C_1 * step_size_array;
    loglog(step_size_array, errors, 'r-o', 'LineWidth', 2)
    hold on;
    loglog(step_size_array,untere_schranke,"g--","Linewidth",1.5);
    grid on;
    ylim([1e-2 1e2]);
    xlabel("logarithmische Schrittweiten");
    ylabel("globaler Fehler (logarithmisch)");
    legend("Exp. Euler Fehler","C_1 tau  ,C_1 ∊ ℝ ","Location","northwest");
