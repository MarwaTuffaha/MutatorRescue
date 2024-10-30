% Code created by Loïc Marrec
% Edited by Marwa Z Tuffaha

function [pr_W,pr_M,pr_all] = Gillespie_fct(Nit, n, m, theta, gW, XW_i, gRW, XRW_i, gM, XM_i, gRM, XRM_i, K, nu, r, e, F, sim)
    
    rng("shuffle")
    mu=r*nu;

    XRW_list = NaN(Nit,1);
    XRM_list = NaN(Nit,1);

    for i = 1 : Nit


        % Initialization

        XW = XW_i;              % Initialization of the number of W individuals (WT)
        XRW = XRW_i;              % Initialization of the number of RW individuals (Rescue from WT)
        XM = XM_i;              % Initialization of the number of M individuals (mutator)
        XRM = XRM_i;              % Initialization of the number of M individuals (Rescue from mutator)
        t = 0;                  % Initialization of time
        cumul = zeros(11, 1);   % To build the sampling tower
        
        mu_M=0; %Mutator emergence rate
        mu_RM=0; %Rescue mutant rate from mutator
        mu_RW=0; %Rescue mutant rate from WT
        if sim==0 % only WT mutates into rescue (no mutator population)
            mu_RW=mu; %Rescue mutant rate from WT
        elseif sim==1 % Mutator population emerges but only WT mutates into rescue
            mu_M=e*mu; %Mutator emergence rate
            mu_RW=mu; %Rescue mutant rate from WT
        elseif sim==2 % Mutator population emerges and only mutator mutates into rescue
            mu_M=e*mu; %Mutator emergence rate
            mu_RM=F*mu; %Rescue mutant rate from mutator
        elseif sim==3 % RWoth mutator and WT mutate
            mu_M=e*mu; %Mutator emergence rate
            mu_RW=mu; %Rescue mutant rate from WT
            mu_RM=F*mu; %Rescue mutant rate from mutator        
        end

        while XW ~= 0 && XRW<K && XRM<K

            fW = max(0,1-nu)*sigm(t, theta, n);
            fRW = max(0,1-nu)*(1-sigm(t, theta, m));
            fM = max(0,1-F*nu)*sigm(t, theta, n);
            fRM = max(0,1-F*nu)*(1-sigm(t, theta, m));            
            
            % Mompute the transition rates

            T_W_rep_no_mut = fW*(1-(XW+XRW+XM+XRM)/K)*(1-mu_RW-mu_M)*XW; 
            T_W_rep_mut_rescue = fW*(1-(XW+XRW+XM+XRM)/K)*mu_RW*XW;
            T_W_rep_mut_mutator = fW*(1-(XW+XRW+XM+XRM)/K)*mu_M*XW;
            T_W_death = gW*XW;
            T_RW_rep = fRW*(1-(XW+XRW+XM+XRM)/K)*XRW;
            T_RW_death = gRW*XRW;
            T_M_rep_no_mut = fM*(1-(XW+XRW+XM+XRM)/K)*(1-mu_RM)*XM; 
            T_M_rep_mut_rescue = fM*(1-(XW+XRW+XM+XRM)/K)*mu_RM*XM;
            T_M_death = gM*XM;
            T_RM_rep = fRM*(1-(XW+XRW+XM+XRM)/K)*XRM;
            T_RM_death = gRM*XRM;            

            T = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep+T_RW_death+T_M_rep_no_mut+T_M_rep_mut_rescue+T_M_death+T_RM_rep+T_RM_death;

            % Mompute tau and then the new time

            r1 = rand;
            tau = 1/T*log(1/r1);
            t = t+tau;

                % Build a sampling tower

                ir2 = 1;
                r2 = rand;
                cumul(1) = T_W_rep_no_mut;
                cumul(2) = T_W_rep_no_mut+T_W_rep_mut_rescue;
                cumul(3) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator;
                cumul(4) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death;
                cumul(5) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep;
                cumul(6) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep+T_RW_death;
                cumul(7) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep+T_RW_death+T_M_rep_no_mut;
                cumul(8) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep+T_RW_death+T_M_rep_no_mut+T_M_rep_mut_rescue;
                cumul(9) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep+T_RW_death+T_M_rep_no_mut+T_M_rep_mut_rescue+T_M_death;
                cumul(10) = T_W_rep_no_mut+T_W_rep_mut_rescue+T_W_rep_mut_mutator+T_W_death+T_RW_rep+T_RW_death+T_M_rep_no_mut+T_M_rep_mut_rescue+T_M_death+T_RM_rep;
                cumul(11) = T;

                % Determine which reaction occurs and compute the number of
                % individuals

                while cumul(ir2) < r2*T

                    ir2 = ir2+1;

                end

                if ir2 == 1

                    XW = XW+1;
                    
                elseif ir2 == 2
                    
                    XRW = XRW+1;

                elseif ir2 == 3
                        
                    XM = XM+1;

                elseif ir2 == 4

                    XW = XW-1;

                elseif ir2 == 5

                    XRW = XRW+1;

                elseif ir2 == 6

                    XRW = XRW-1;

                elseif ir2 == 7

                    XM = XM+1;       

                elseif ir2 == 8

                    XRM = XRM+1;                      

                elseif ir2 == 9

                    XM = XM-1; 

                elseif ir2 == 10

                    XRM = XRM+1;

                elseif ir2 == 11

                    XRM = XRM-1;                      

                end
   

            
        end

        XRW_list(i) = XRW;
        XRM_list(i) = XRM;
        
    end
    
    XRWM_list=XRW_list+XRM_list;

    % Compute the rescue probability

    pr_W = length(XRW_list(XRW_list ~= 0))/Nit; % by WT
    pr_M = length(XRM_list(XRM_list ~= 0))/Nit; % by mutator
    pr_all = length(XRWM_list(XRWM_list ~= 0))/Nit; % by all
   
end
