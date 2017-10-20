function [CoolingProfile] = BiotCoolingProfile(z,t,Biot,Kappa)
% This function calculates the Biot number-dependent cooling profile for a
% given set of conditions.
%
% Inputs: z     - The depths at which the profile is calculated.
%         t     - The times at which the profile is calculated.
%         Biot  - The Biot number, with Biot = inf representing a perfectly conducting boundary.
%         Kappa - Thermal/Saline Diffusivity.
%
% (09/11/15)

    CoolingProfile = erf(z/sqrt(4*Kappa.*t)) + ...
        exp(-z.^2./(4*Kappa.*t)).*erfcx(z/sqrt(4*Kappa.*t) + Biot.*sqrt(Kappa.*t));
end

