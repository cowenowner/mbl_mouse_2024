function [ alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] = PR_dend_gating( VmD, Ca )
%PR_dend_gating returns the rate constants for the dendritic gating variables
%of the Pinsky-Rinzel model, as a function of the membrane potential and
%dendritic calcium concentration.
%
%   [ alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] = PR_dend_gating( VmD, Ca )
%   Input VmD should be dendritic membrane potential (a scalar or a vector of
%   values).
%   Input Ca should be dendritic calcium concentration (a scalar or a
%   vector of values).
%   Returned voltage-dependent rate constant arrays
%   (alpha_mca, beta_mca, alpha_kca, beta_kca)
%   are each of the same size
%   as the input array of membrane potentials.
%
%   Returned calcium-dependent rate constants (alpha_kahp and beta_kahp)
%   are each the same size as the input array of calcium values.
%
%   alpha_mca is calcium activation rate constant
%   beta_mca is calcium deactivation rate constant
%   alpha_kca is calcium-dependent potassium activation rate constant
%   beta_kca is calcium-dependent potassium deactivation rate constant
%   alpha_kahp is after-hyperpolarization activation rate constant
%   beta_kahp is after-hyperpolarization deactivation rate constant


alpha_mca = 1600./( 1+exp(-72*(VmD-0.005)) );

beta_mca = ( VmD == -0.0089 ).*20/0.2 ...
    + ( VmD ~= -0.0089 ) ...
    .*20e3.*(VmD+0.0089)./(exp(200*(VmD+0.0089))-1);


alpha_kca = 2e3*exp(-(0.0535+VmD)/0.027).*(VmD>-0.010) ...
    + exp( (VmD+0.050)/0.011 -(VmD+0.0535)/0.027 )/0.018975 ...
    .* (VmD <= -0.010);
beta_kca = (2e3*exp(-(0.0535+VmD)/0.027)-alpha_kca) ...
    .* (VmD <= -0.010);


alpha_kahp = min(20,20e3*Ca);

beta_kahp = 4*ones(size(alpha_kahp));

end

