% To calculate the vector potential in a 3D vox following Yang et. al, 2013

%We use a fourth-order accurate fast direct scheme according to Boisvert (1984), 
%which is included, e.g., in the IMSL Numerical Libraries, to solve the Poisson equation
%(Equation (8)).

% This requires the partial differential equation toolbox
model = createpde();