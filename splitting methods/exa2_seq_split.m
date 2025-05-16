%% figure 6.2

%%  sandbox of variables
    A = [-2 0; 0 3];
    B = [0 2; -4 0];
    u0 = [2; 1];
    step_size_array = 2.^(0:7)*0.001;
    T = 2;

%% pre-allocating
    fehler = zeros(1, length(step_size_array));

%% step_size iteration 
    k =1 ;
    for step_size = step_size_array

        N = round(T/step_size);
        expB = expm(step_size * B);
        expA = expm(step_size * A);
        u_approx = u0;
        fehler_pro_t = zeros(1,N);
        
        %% sequential splitting
        for p = 1:N
    
            first_split  = expA * u_approx;
            second_split = expB * first_split;
            u_approx = second_split ;
            
        %% exact solution
            u_exact = expm(p*step_size* (A+B))*u0;
    
            fehler_pro_t(p) = norm(u_approx - u_exact,inf);

        end
        
    %% global error
        fehler(k) = norm(fehler_pro_t,inf);
        k = k+1;
    end

%% Plotting 
    figure("Name","error plot sequential splitting method for matrix", ...
        "NumberTitle","off"); 
    C_1 = 50;      
    untere_schranke    = C_1 * (step_size_array);
    loglog(step_size_array,fehler,"r-o","Linewidth",2);
    hold on;
    loglog(step_size_array,untere_schranke,"g--","Linewidth",1.5);
    grid on;
    ylim([1e-2 1e2]);
    xlabel("Schrittweiten (\tau, logarithmisch");
    ylabel("globaler Fehler (logarithmisch)");
    legend("Seq. Splitting Fehler","C_1 tau  ,C_1 ∊ ℝ ","Location","northwest");