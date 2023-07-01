% Project1 Computational Fluid Dynamics
% Profossor : DR.Naderan 
% Student : Vahid Eftekhari Khorasani
% Student Number: 401126104
% 2D Lid-driven cavity solved with SIMPLE algoritm 

clc;
clear ;
close all;

Reynolds=input(' Enter Reynold_number: 1-1: 2-10: 3-100: 4-500 ');
[Re]= Reynolds_Number(Reynolds);

if Re < 100
        N = [33 41 49 57 65 ];
    else
        N=[129 137 145 153 161 ];
end

Grid   = length(N);
Fact   = zeros(4,Grid);
U_Mean = zeros(4,Grid);
H      = zeros(1,Grid);

%% Describe the domain of the problem
for ii=1:Grid    
    N_nods = N(ii);  % Number of nods
Wall_length = 1;
    h_grid = Wall_length/(N_nods-1); 
        x = 0:h_grid:Wall_length; %X domain span
    y = 0:h_grid:Wall_length; %Y domain span 
    
    nu = 1/Re;
    Error1=zeros(N_nods-1,N_nods-1,5);

    % Relaxation factors for Convergence
    % Under-relaxation factor (alpha)
    alpha = 0.8;    % for velocity
    alpha_p = 0.8;  % for pressure
    U_top = 1;
    %% Initializing the variables
    % End_collocated_variables
    U_end = zeros(N_nods,N_nods);
    V_end = zeros(N_nods,N_nods);
    P_end = zeros(N_nods,N_nods);
    U_end(1,:) = U_top;
    
    %Staggered variables 
    U       = zeros(N_nods+1,N_nods);
    U_stars = zeros(N_nods+1,N_nods);
    d_e     = zeros(N_nods+1,N_nods);
    V       = zeros(N_nods,N_nods+1);
    V_stars = zeros(N_nods,N_nods+1);
    d_n     = zeros(N_nods,N_nods+1);
    P       = zeros(N_nods+1,N_nods+1);
    P_stars = zeros(N_nods+1,N_nods+1);
    
    P(N_nods+1,N_nods+1) = 1;
    P_stars(N_nods+1,N_nods+1) = 1;
    p_c = zeros(N_nods+1,N_nods+1);
    b = zeros(N_nods+1,N_nods+1);
    U(1,:) = 2*U_top;
    
    U_new = zeros(N_nods+1,N_nods);
    V_new = zeros(N_nods,N_nods+1);
    P_new = zeros(N_nods+1,N_nods+1);
    P_new(N_nods+1,N_nods+1)=1;
    U_new(1,:) = 2*U_top;
    Data_table = zeros(17,2);
    %% Solving the governing equations
    error = 1;
    iterations = 0;
    error_con = 1e-7; %final required error residual
    figure(1);        %for error monitoring
    
    while error > error_con
        % X direction momentum equation inside volume
        for i = 2:N_nods
            for j = 2:N_nods - 1
                u_E = 0.5*(U(i,j) + U(i,j+1));
                u_W = 0.5*(U(i,j) + U(i,j-1));
                v_N = 0.5*(V(i-1,j) + V(i-1,j+1));
                v_S = 0.5*(V(i,j) + V(i,j+1));
                
                a_E = -0.5*u_E*h_grid + nu;
                a_W = 0.5*u_W*h_grid + nu;
                a_N = -0.5*v_N*h_grid + nu;
                a_S = 0.5*v_S*h_grid + nu;
                
                a_e = 0.5*u_E*h_grid - 0.5*u_W*h_grid + 0.5*v_N*h_grid - 0.5*v_S*h_grid + 4*nu;
                
                A_e = -h_grid;
                d_e(i,j) = A_e/a_e;
                
                U_stars(i,j) = (a_E*U(i,j+1) + a_W*U(i,j-1) + a_N*U(i-1,j) + a_S*U(i+1,j))/a_e + d_e(i,j)*(P(i,j+1) - P(i,j));
            end
         end
    
    % X direction momentum equation Boundary Condition
    U_stars(1,:) = 2.*U_top - U_stars(2,:);
    U_stars(N_nods + 1,:) = -U_stars(N_nods,:);
    U_stars(2:N_nods,1) = 0;
    U_stars(2:N_nods,N_nods) = 0;
    
    % Y direction momentum equation inside volume
    for i = 2:N_nods - 1
        for j = 2:N_nods
            u_E = 0.5*(U(i,j) + U(i+1,j));
            u_W = 0.5*(U(i,j-1) + U(i+1,j-1));
            v_N = 0.5*(V(i-1,j) + V(i,j));
            v_S = 0.5*(V(i,j) + V(i+1,j));
            
            a_E = -0.5*u_E*h_grid + nu;
            a_W = 0.5*u_W*h_grid + nu;
            a_N = -0.5*v_N*h_grid + nu;
            a_S = 0.5*v_S*h_grid + nu;
            
            a_n = 0.5*u_E*h_grid - 0.5*u_W*h_grid + 0.5*v_N*h_grid - 0.5*v_S*h_grid + 4*nu;
            
            A_n = -h_grid;
            d_n(i,j) = A_n/a_n;
            
            V_stars(i,j) = (a_E*V(i,j+1) + a_W*V(i,j-1) + a_N*V(i-1,j) + a_S*V(i+1,j))/a_n + d_n(i,j)*(P(i,j) - P(i+1,j));
        end
    end
    
    % Y direction momentum equation Boundary Condition
    V_stars(:,1) = -V_stars(:,2);
    V_stars(:,N_nods + 1) = -V_stars(:,N_nods);
    V_stars(1,2:N_nods) = 0;
    V_stars(N_nods,2:N_nods) = 0;
    
    %  Modifying p_c to zero (correction)
    p_c(1:N_nods+1,1:N_nods+1)=0;
    
    % Continuity equation __ pressure correction inside volume
    for i = 2:N_nods
        for j = 2:N_nods
            a_E = -d_e(i,j)*h_grid;
            a_W = -d_e(i,j-1)*h_grid;
            a_N = -d_n(i-1,j)*h_grid;
            a_S = -d_n(i,j)*h_grid;
            a_P = a_E + a_W + a_N + a_S;
            b(i,j) = -(U_stars(i,j) - U_stars(i,j-1))*h_grid + (V_stars(i,j) - V_stars(i-1,j))*h_grid;
            
            p_c(i,j) = (a_E*p_c(i,j+1) + a_W*p_c(i,j-1) + a_N*p_c(i-1,j) + a_S*p_c(i+1,j) + b(i,j))/a_P;
        end
    end
    
    % Modifying the pressure field
    for i = 2:N_nods
        for j = 2:N_nods
            P_new(i,j) = P(i,j) + alpha_p*p_c(i,j);
        end
    end
    
    % Continuity equation Boundary Condition
    P_new(1,:) = P_new(2,:);
    P_new(N_nods + 1,:) = P_new(N_nods,:);
    P_new(:,1) = P_new(:,2);
    P_new(:,N_nods + 1) = P_new(:,N_nods);
    
    % Modifying the velocity
    for i = 2:N_nods
        for j = 2:N_nods - 1
            U_new(i,j) = U_stars(i,j) + alpha*d_e(i,j)*(p_c(i,j+1) - p_c(i,j));
        end
    end
    
    % X direction momentum equation Boundary
    U_new(1,:) = 2*U_top - U_new(2,:);
    U_new(N_nods + 1,:) = -U_new(N_nods,:);
    U_new(2:N_nods,1) = 0;
    U_new(2:N_nods,N_nods) = 0;
    
    for i = 2:N_nods - 1
        for j = 2:N_nods
            V_new(i,j) = V_stars(i,j) + alpha*d_n(i,j)*(p_c(i,j) - p_c(i+1,j));
        end
    end
    
    % Y direction momentum equation Boundary Condition
    V_new(:,1) = -V_new(:,2);
    V_new(:,N_nods + 1) = -V_new(:,N_nods);
    V_new(1,2:N_nods) = 0;
    V_new(N_nods,2:N_nods) = 0;
            
    
    % Obtain Residual for error measure
    error = 0;
    for i = 2:N_nods
        for j = 2:N_nods
            error = error + abs(b(i,j));
        end
    end
    
    % variables in the last of iteration
    U = U_new;
    V = V_new;
    P = P_new;
    iterations = iterations + 1;
    Error1(iterations,ii)=error;
    end

    % transform the staggered variables to collocated variables , after
    % convergence
    for i = 1:N_nods
        for j = 1:N_nods
            U_end(i,j) = 0.5*(U(i,j) + U(i+1,j));
            V_end(i,j) = 0.5*(V(i,j) + V(i,j+1));
            P_end(i,j) = 0.25*(P(i,j) + P(i,j+1) + P(i+1,j) + P(i+1,j+1));
        end
    end


%% Plot Setting and details
x_int = ((1:N_nods)-1).*h_grid;
y_int = 1-((1:N_nods)-1).*h_grid;
[X,Y] = meshgrid(x_int,y_int);

%plot U velocity_X
figure(1);
subplot(2,3,ii)
contourf(X,Y,U_end, 21, 'LineStyle', 'none')

colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title(sprintf(' U velocity X direction (m/s) for N= %d ' ,round(N(ii))))

% plot U velocity_X
figure(2);
subplot(2,3,ii)
contourf(X,Y,V_end, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title(sprintf(' V velocity in Y direction (m/s) for N= %d ' ,round(N(ii))))

%Pressure(pa)
figure(3);
subplot(2,3,ii)
contourf(X,Y,P_end, 49, 'LineStyle' , "none")
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title(sprintf('Pressure(pa) for N= %d ' ,round(N(ii))))

% plot Total velocity vector
figure(4);
subplot(2,3,ii)
hold on
grid on
ploth = quiver(X, Y, U_end, V_end, 5, 'k');
shading interp;
title(sprintf('Total velocity vector for N= %d ',round(N(ii))))
axis equal

% Error 
ite=1:iterations;

figure(5);
subplot(2,3,ii)
plot(ite,Error1(:,ii))
hold on
xlabel('Iterations')
ylabel('Residual Error')
title(sprintf('Residual Error(Iterations) for N= %d ' ,round(N(ii))))
%% Validation with Dr. Qia's article for RE == 100
if Re == 100
    figure(6);
    plot(U_end(:,(N_nods+1)/2), 1-y, 'LineWidth', 1)
    Data_table = xlsread('C:\Users\LENOVO\OneDrive\Desktop\CFD\Project CFD\Lid-drive cavity flow\Ghia_Re_100.csv','A2:B18');
    y_ghia = Data_table(:,1);
    u_ghia = Data_table(:,2);
    
    figure(6); hold on
    plot(u_ghia, y_ghia, 'o', 'LineWidth', 1)
    title(sprintf('Residual Error(Iterations) for N= %d ' ,round(N(ii))))
    xlabel('u')
    ylabel('y')
    legend('Numerical', 'Benchmark', 'location', 'southeast')
end


    H(ii)=h_grid;
    Fact(:,ii)=[find(y==0.25);find(y==0.5);find(y==0.625);find(y==0.75)];
for i=1:4
    U_Mean(i,ii)=U_end((N_nods+1)/2,Fact(i,ii));
end

end 

Eror =Error(U_Mean,Grid);
Delta_y = Error_Slope( Eror,H,Grid );

%% Plot the relative error value in terms of H
COL=['s','*','d','o'];

for i=1:4
    figure(7)
    subplot(1,2,1)
    loglog(H(1:end-1),abs(Eror(i,:)),sprintf(COL(i)),'linewidth',1.5);
    loglog(H(1:end-1),abs(Eror(i,:)),'linewidth',1.5);
    polyfit(H(1:end-1), abs(Eror(i,:)), 1);
    title('Solution convergence ');
    ylabel('successive error');
    xlabel('h');
    hold on
    grid on
    
    %% PLot the error slope in terms of H
    subplot(1,2,2)
    loglog(H(1:end-2),abs(Delta_y(i,:)),'linewidth', 2);
    title('Error slope ');
    ylabel('Slope');
    xlabel('h');
    hold on
    grid on
end
legend('y=0.25','y=0.5','y=0.75','y=0.875')
