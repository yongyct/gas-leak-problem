function [] = gas_leak_problem()

%to calculate for 5 species
for j = 1:5                 

%initial guess of friction factor
f(1,j) = 0.02;             

%defining given constants
L = 8;
D = 3e-2;
gamma = [1.4; 1.316; 1.2461; 1.32; 1.26];
R = 8.314;
T = 293;
P1 = 3e5;
M = [29e-3; 17.031e-3; 28.05e-3; 16.04e-3; 64.066e-3];
c = sqrt(gamma.*R.*T./M);
rho = P1.*M./R/T;
mu = [1.82e-5; 0.99e-5; 1.03e-5; 1.1e-5; 1.26e-5];
rr = 0.00005;

%setting options for solvers
options = optimset('TolX', 1e-20);

%assuming choked conditions i.e. Ma2 = 1
Ma2(j) = 1;
    
%iteration for convergence of friction factor
for i = 1:10
    sma1 = @solnma1;
    ma1(i,j) = fzero(sma1, 0.5, options, f(i,j), L , D, gamma(j));
    v(i,j) = ma1(i,j)*c(j);
    re1(i,j) = rho(j)*v(i,j)*D/mu(j);
    f(i+1,j) = (-1.8*log10((rr/3.7)^1.11 + 6.9/re1(i,j)))^(-2);
    if(abs(f(i+1,j)-f(i,j)) < 0.0001)
        break
    end
end

%calculating answers for question
Ma1(j) = ma1(end,j);
P2(j) = 3e5*ma1(end,j)...
    *sqrt((1+(gamma(j)-1)/2*(ma1(end,j))^2)/(1+(gamma(j)-1)/2));
T2(j) = T*(1+(gamma(j)-1)/2*(ma1(end,j))^2)/(1+(gamma(j)-1)/2);
G(j) = rho(j)*v(end,j);

%for unchoked conditions
if(P2(j)<1e5)
    ma1u(1,j) = 0.3;
    a(1,j) = 1e-7;
    
%iteration for convergence of Ma1
    for i = 1:10
        vu(i,j) = ma1u(i,j)*c(j);
        re1u(i,j) = rho(j)*vu(i,j)*D/mu(j);
        fu(i,j) = (-1.8*log10((rr/3.7)^1.11 + 6.9/re1u(i,j)))^(-2);
        sma2 = @solnma2;
        ma2u(i,j) = ...
            fzero(sma2, 1, options, ma1u(i,j), fu(i,j), L, D, gamma(j));
        P2u(i,j) = P1*ma1u(i,j)/ma2u(i,j)*...
            sqrt((1+(gamma(j)-1)/2*(ma1u(i,j))^2)/...
            (1+(gamma(j)-1)/2*(ma2u(i,j))^2));
        if (abs(P2u(i,j) - 1e5) <= 1000)
            break            
        end
        if (i == 1)
        else
            a(i,j) = 0.29935*(i-1)*a(i-1,j);
        end                    
            ma1u(i+1,j) = ma1u(i,j) + a(i,j)*(P2u(i,j)-1e5);   
    end
    T2u(j) = T*(1+(gamma(j)-1)/2*(ma1u(end,j))^2)/...
        (1+(gamma(j)-1)/2*(ma2u(end,j))^2); 
    Gu(j) = rho(j)*vu(end,j);
    
%finalising answers
    Ma1(j)= ma1u(end,j);
    Ma2(j) = ma2u(end,j);
    T2(j) = T2u(j);
    P2(j) = P2u(end,j);
    G(j) = Gu(j);
    
end

end

%%
%function relating Ma1 to f
function solver = solnma1(ma1,f,L,D,gamma)
solver = (gamma+1)/2*log((1+(gamma-1)/2*ma1^2)/(ma1^2*(1+(gamma-1)/2)))...
    - 1/(ma1^2) + 1 + gamma*f*L/D;
end

%%
%function relating Ma1, Ma2 and f
function solver = solnma2(ma2,ma1,f,L,D,gamma)
solver = (gamma+1)/2*...
    log((ma2^2*(1+(gamma-1)/2*ma1^2))/(ma1^2*(1+(gamma-1)/2*(ma2^2))))...
    - 1/(ma1^2) + 1/(ma2^2) + gamma*f*L/D;
end

%%
%Tabulating answers
Rowname ={'Air';'Ammonia';'Ethylene';'Methane';'Sulfur Dioxide'};
Ma_In = Ma1';Ma_Out = Ma2';T_Out = T2';P_Out = P2'; Gas_Flow = G';
TBL = table(Ma_In,Ma_Out,T_Out,P_Out,Gas_Flow,'RowNames',Rowname);
TBL1 = table(Ma_In,Ma_Out,'RowNames',Rowname)
TBL2 = table(T_Out,P_Out,Gas_Flow,'RowNames',Rowname)

%%
%Appendix

for j = 1:5
    
    switch j
        case 1
            Air_Choked = ...
                table(f(1:3,j), ma1(:,j), re1(:,j), f(2:4,j));
            Air_Choked.Properties.Description = 'Air (Choked)';
            Air_Choked.Properties.VariableNames = ...
                {'f_old' 'Ma_1' 'Re_1' 'f_new'};
            Air_Unchoked = ...
                table(ma1u(:,j), re1u(:,j), fu(:,j), ma2u(:,j), P2u(:,j));
            Air_Unchoked.Properties.Description = 'Air (Unchoked)';
            Air_Unchoked.Properties.VariableNames = ...
                {'Ma_1' 'Re_1' 'f' 'Ma_2' 'P_2'};           
        case 2
            Ammonia_Choked = table(f(1:3,j), ma1(:,j), re1(:,j), f(2:4,j));
            Ammonia_Choked.Properties.Description = 'Ammonia (Choked)';
            Ammonia_Choked.Properties.VariableNames = ...
                {'f_old' 'Ma_1' 'Re_1' 'f_new'};
        case 3
            Ethylene_Choked=table(f(1:3,j), ma1(:,j), re1(:,j), f(2:4,j));
            Ethylene_Choked.Properties.Description = 'Ethylene (Choked)';
            Ethylene_Choked.Properties.VariableNames = ...
                {'f_old' 'Ma_1' 'Re_1' 'f_new'};
        case 4
            Methane_Choked = table(f(1:3,j), ma1(:,j), re1(:,j), f(2:4,j));
            Methane_Choked.Properties.Description = 'Methane (Choked)';
            Methane_Choked.Properties.VariableNames = ...
                {'f_old' 'Ma_1' 'Re_1' 'f_new'};
        case 5
            Sulfur_Dioxide_Choked =...
                table(f(1:3,j), ma1(:,j), re1(:,j), f(2:4,j));
            Sulfur_Dioxide_Choked.Properties.Description =...
                'Sulfur Dioxide (Choked)';
            Sulfur_Dioxide_Choked.Properties.VariableNames = ...
                {'f_old' 'Ma_1' 'Re_1' 'f_new'};
    end  
end

Air_Choked
Air_Unchoked
Ammonia_Choked
Ethylene_Choked
Methane_Choked
Sulfur_Dioxide_Choked

end