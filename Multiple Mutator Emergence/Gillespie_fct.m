% Code created by Loïc Marrec
% Edited by Marwa Z Tuffaha

function [pr_W,pr_M,pr_all] = Gillespie_fct(Nit, n, m, theta, gW, XW_i, gRW, XRW_i, gM, XM_i, gRM, XRM_i, K, nu, r, e, F_vals, sim)
    
    rng("shuffle")
    mu=r*nu;
    F_N=length(F_vals);
    mu_M_i= e/F_N *mu; % emergence rate for each mutator

    XRW_list = NaN(Nit,1);
    XRM_list = NaN(Nit,1);

    %i_counter_print=500;
    
    for i = 1 : Nit

        % Initialization

        XW = XW_i;              % Initialization of the number of WT individuals
        XRW = XRW_i;              % Initialization of the number of Rescue from WT individuals
        XM = XM_i*ones(1,F_N);              % Initialization of the number of mutators individuals
        XRM = XRM_i*ones(1,F_N);             % Initialization of the number of rescue from mutators individuals
        t = 0;                  % Initialization of time
        
        mu_M=zeros(1,F_N); %Mutator emergence rate
        mu_RM= zeros(1,length(F_vals)); %Rescue mutant rate from mutator
        mu_RW=0; %Rescue mutant rate from WT
        if sim==0 % only WT mutates into rescue (no mutator population)
            mu_RW=mu; %Rescue mutant rate from WT
        elseif sim==1 % Mutator population emerges but only WT mutates into rescue
            mu_M=mu_M_i*ones(1,length(F_vals)); %Mutator emergence rate
            mu_RW=mu; %Rescue mutant rate from WT
        elseif sim==2 % Mutator population emerges and only mutator mutates into rescue
            mu_M=mu_M_i*ones(1,length(F_vals)); %Mutator emergence rate
            mu_RM=F_vals.*ones(1,length(F_vals))*mu; %Rescue mutant rate from mutator
        elseif sim==3 % Both mutator and WT mutate
            mu_M=mu_M_i*ones(1,length(F_vals)); %Mutator emergence rate
            mu_RW=mu; %Rescue mutant rate from WT
            mu_RM=F_vals.*ones(1,length(F_vals))*mu; %Rescue mutant rate from mutator        
        end

        while XW ~= 0 && XRW<K && sum(XRM)<K

            fW = (1-nu)*sigm(t, theta, n);
            fRW = (1-nu)*(1-sigm(t, theta, m));
            fM = max(0,1-F_vals.*nu).*sigm(t, theta, n);
            fRM = max(0,1-F_vals.*nu).*(1-sigm(t, theta, m));            
            
            % Compute the transition rates

            log_fac=(1-(XW+XRW+sum(XM)+sum(XRM))/K);

            T_W_rep_no_mut = fW*log_fac*(1-mu_RW-sum(mu_M))*XW; 
            T_W_rep_mut_rescue = fW*log_fac*mu_RW*XW;
            T_W_rep_mut_mutator = fW*log_fac.*mu_M*XW;
            T_W_death = gW*XW;
            T_RW_rep = fRW*log_fac*XRW;
            T_RW_death = gRW*XRW;
            T_M_rep_no_mut = fM.*log_fac.*(1-mu_RM).*XM; 
            T_M_rep_mut_rescue = fM.*log_fac.*mu_RM.*XM;
            T_M_death = gM*XM;
            T_RM_rep = fRM.*log_fac.*XRM;
            T_RM_death = gRM*XRM;            

            T = T_W_rep_no_mut+T_W_rep_mut_rescue+sum(T_W_rep_mut_mutator)+T_W_death+T_RW_rep+T_RW_death+sum(T_M_rep_no_mut)+sum(T_M_rep_mut_rescue)+sum(T_M_death)+sum(T_RM_rep)+sum(T_RM_death);
            
            % sampling tower
            Rates = [T_W_rep_no_mut,T_W_rep_mut_rescue,T_W_death,T_RW_rep,T_RW_death,T_W_rep_mut_mutator,T_M_rep_no_mut,T_M_rep_mut_rescue,T_M_death,T_RM_rep,T_RM_death];
            cumul=cumsum(Rates);
            N_events=length(cumul);

            % Compute tau and then the new time

            r1 = rand;
            tau = 1/T*log(1/r1);
            t = t+tau;

                % Build a sampling tower

                ir2 = 1;
                r2 = rand;

                % Determine which reaction occurs and compute the number of individuals

                while cumul(ir2) < r2*T
                    ir2 = ir2+1;
                    if ir2>N_events
                        disp(">N_events")
                    end

                end

                if ir2 == 1 % WT replicates no mutation

                    XW = XW+1;
                    
                elseif ir2 == 2 % WT replicates and mutates to rescue
                    
                    XRW = XRW+1;

                elseif ir2 == 3 % WT death

                    XW = XW-1;

                elseif ir2 == 4 % Rescue from WT replicates

                    XRW = XRW+1;

                elseif ir2 == 5 % Rescue from WT death

                    XRW = XRW-1;

                elseif ir2 <= 5+F_N % WT replicates and mutates to mutator
                    
                    j=mod(ir2-5,F_N+1);
                        XM(j) = XM(j)+1;


                elseif ir2 <=5+2*F_N % Mutator replicates no mutation
                    j=mod(ir2-5-F_N,F_N+1);
                    XM(j) = XM(j)+1;     

                elseif ir2 <= 5+3*F_N % mutator replicates and mutates into rescue
                    j=mod(ir2-5-2*F_N,F_N+1);
                    XRM(j) = XRM(j)+1;                     

                elseif ir2 <=5+4*F_N % Mutator death
                    j=mod(ir2-5-3*F_N,F_N+1);
                    XM(j) = XM(j)-1;   

                elseif ir2 <=5+5*F_N % Rescue from mutator replicates
                    j=mod(ir2-5-4*F_N,F_N+1);
                    XRM(j) = XRM(j)+1; 

                elseif ir2 <=5+6*F_N % Rescue from mutator death
                    j=mod(ir2-5-5*F_N,F_N+1);
                    XRM(j) = XRM(j)-1;                     

                end
            
        end

        XRW_list(i) = XRW;
        XRM_list(i) = sum(XRM);
        
    end
    
    XR_list=XRW_list+XRM_list;

    % Count how many times it was successfully rescued
     
    pr_W = length(XRW_list(XRW_list ~= 0))/Nit;
    pr_M = length(XRM_list(XRM_list ~= 0))/Nit;
    pr_all= length(XR_list(XR_list ~= 0))/Nit;
   
end
