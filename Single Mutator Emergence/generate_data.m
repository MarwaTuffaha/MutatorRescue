% Code created by Loïc Marrec
% Edited by Marwa Z Tuffaha

% Parameters
n = 2;                          % Hill coefficient for Wildtype otr mutator
m = n;                          % Hill coefficient for rescue mutants
theta = 500;                    % Inflection point
gW = 0.1;                       % Death rate of W microbes
XW_i = 9000;                    % Initial number of W microbes
gRW = 0.1;                      % Death rate of rescue mutants originating from W microbes
XRW_i = 0;                      % Initial number of rescue mutants originating from W microbes
gM = 0.1;                       % Death rate of mutator microbes
XM_i = 0;                       % Initial number of mutator microbes
gRM = 0.1;                      % Death rate of rescue mutants originating from mutator microbes
XRM_i = 0;                      % Initial number of rescue mutants originating from mutator microbes
K = 1e4;                        % Carrying capacity
Nit = 500;                      % Number of stochastic realizations 
nu=0.001;                       % Wildtype mutation rate
F=5;                            % Mutator strength

r=2*10^-4;                      % Fraction of rescue mutations out of all available mutations
                                % Wildtype's mutation rate at the rescue allele is mu=r*nu

e=1000;                         % Mutator emergence factor
                                % Mutator emergence rate is mu_M=e*mu

                                                            
% parameter deciding what kind of simulation we want
% 0=only WT mutates into rescue (no mutator population)
% 1=Mutator population emerges but only WT mutates into rescue
% 2=Mutator population emerges and only mutator mutates into rescue
% 3= both mutate
sim=3;


[pr_W,pr_M,pr_all] = Gillespie_fct(Nit, n, m, theta, gW, XW_i, gRW, XRW_i, gM, XM_i, gRM, XRM_i, K, nu, r, e, F, sim);




