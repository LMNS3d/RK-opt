function [stability] = internal_stability_butcher(A,b,c,spectrum,one_step_dt)
%function [stability] = internal_stability(rk_coeffs_file,rk_form,spectrum_file)
%
%
% This function check if the Runge-Kutta scheme given in the Butcher form
% is internally stable.
%
%%%%%%%%%%%%%%%%%%%%%%%%


% Check intermediate stability functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number os stages
s = length(A);

% Combine A with b^T to get the characterizing matrix
K = [A;b'];

% Construct the intermediate stability functions \psi_j (where j is the 
% index of the stage), verify that the scaled spectrum (spectrum * one_step_dt = z) 
% lies inside it and compute the linearly stable time step for such a 
% intermediate stability function  
for i = 1:s
    % Get generalized A and b^T
    [gA,gbT] = sub(K,i)
    
    % Cosntruct identity matrix of the right size
    Id =  eye(length(gA));
    
    % Construct unit vector of the right size
    ge = ones(length(gA),1);
    
    % Define symbolic stability function
    syms z
    stab_function = matlabFunction(1 + gbT*z*(Id - gA*z)^(-1)*ge);
    
    [delta_t_inter] = stable_time_step(stab_function,1);


    % Plot absolute stability region in the complex plane
    stabPlot = figure;
    hand = gca;
    x = linspace(-3,1,200);
    y = linspace(3,-3,200);
    [X,Y] = meshgrid(x,y);
    Z = X + Y*sqrt(-1);
    
    stab_function_eval = subs(stab_function,{z},{Z});
    
    [c,h] = contour(X,Y,stab_function_eval,[1,1]);
    set(h,'edgecolor','r');

end

end


% =========================================================================


% =========================================================================


function [delta_t] = stable_time_step(stab_fun,spec)
%function 
%
% Bisection method to find the stable time step for a given stability
% function
%%%%%%%%%%%%%%%%%%%%%%%%

% Bisection parameters
eps = 1.e-6;
toll = 1.e-10;

delta_t_min = 0.0; 
delta_t_max = 20.0;
max_z = 0.0;

% Symbolic independent variable
syms z

while (delta_t_max - delta_t_min) > eps
    delta_t = (delta_t_min + delta_t_max)/2.0;
    
    scaled_spec = delta_t*spec;
    
    stab_fun_eval = subs(stab_fun,{z},{scaled_spec});
    
    max_Z = max(abs(stab_fun_eval));
    
    if max_z > (1 + toll)
        delta_t_max = delta_t;
    else
        delta_t_min = delta_t;
    end
end

end
% =========================================================================


% =========================================================================

function [C,w] = sub(K,i)
%function mat = sub(K,i)
%
% Extract from the matrix K, the leading principal submatrix of order i
% K_i = K(1:i,1:i), and the vector of the weights.
%%%%%%%%%%%%%%%%%%%%%%%%

C = K(1:i, 1:i);
w = K(i+1,1:i);

% =========================================================================

end
