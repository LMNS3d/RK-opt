function [stability] = internal_stability(rk_coeffs_file,rk_form,spectrum_file)
%function [stability] = internal_stability(rk_coeffs_file,rk_form,spectrum_file)
%
%
% This function check if the Runge-Kutta scheme uniquely defined by the
% coefficients provided in the rk_coeffs_file is internally stable.
%
% The Runge-Kutta can be specified in the Butcher form or in the modfied 
% Shu-Osher.
% ================================================================

% Load Runge-Kutta coefficients
if strcmp(rk_form,'butcher')
    [rk.A,rk.b,rk.c] = read_butcher_coeffs(rk_coeffs_file);
elseif strcmp(objective,'acc')
    [rk.v,rk.alpah,rk.beta] = read_mSO_coeffs(rk_coeffs_file);
else
    error('Unrecognized Runge-Kutta form');
end

% Load matrix spectrum of ODE system, i.e. du/dt = L*u 
spectrum = read_spectrum(spectrum_file);







end

% =========================================================================


% =========================================================================

function [A,b,c] = read_butcher_coeffs(file)
%function [A,b,c] = read_butcher_coeffs(file)
%


% Load Runge-Kutta coefficients
if strcmp(rk_form,'butcher')
    [rk.A,rk.b,rk.c] = read_butcher_coeffs(rk_coeffs_file);
elseif strcmp(objective,'acc')
    [rk.v,rk.alpah,rk.beta] = read_mSO_coeffs(rk_coeffs_file);
else
    error('Unrecognized Runge-Kutta form');
end

% Load matrix spectrum of ODE system, i.e. du/dt = L*u 
spectrum = read_spectrum(spectrum_file);







end

% =========================================================================