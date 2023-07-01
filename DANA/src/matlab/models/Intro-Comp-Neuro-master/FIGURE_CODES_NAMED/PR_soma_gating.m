function [ alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] = PR_soma_gating( Vm )
%PR_soma_gating returns the rate constants for the somatic gating variables
%of the Pinsky-Rinzel model, as a function of the membrane potential.
%
%   [ alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] = PR_soma_gating( Vm )
%   Input, Vm, should be somatic membrane potential (a scalar or a vector of
%   values).
%   Returned voltage-dependent rate constants are each of the same size as
%   the input membrane potential.
%   alpha_m is sodium activation rate constant
%   beta_m is sodium deactivation rate constant
%   alpha_h is sodium inactivation rate constant
%   beta_h is sodium deinactivation rate constant
%   alpha_n is potassium activation rate constant
%   beta_n is potassium deactivation rate constant


    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. 
    
    alpha_m = ( Vm == -0.0469 ).*320/0.25 ...
        + ( Vm ~= -0.0469 ) ...
        .*( 320*1e3*(Vm+0.0469) )./(1-exp(-250*(Vm+0.0469))) ;
    
    
    beta_m = ( Vm == -19.9 ) .* 280/0.2 ...
        + ( Vm ~= -19.9 ) .* ...
        ( 280e3*(Vm+0.0199) )./(exp(200*(Vm+0.0199))-1);
    

    alpha_h = 128*exp(-(Vm+0.043)/0.018);
    beta_h = 4e3./(1+exp(-200*(Vm+0.020)));
    
    alpha_n = ( Vm == -0.0249 ) .* 16/0.2 ...
        + ( Vm ~= -0.0249 ) .* ...
        ( 16e3*(Vm+0.0249) )./(1-exp(-200*(Vm+0.0249)));
    
    beta_n = 250*exp(-25*(Vm+0.040));    

end

