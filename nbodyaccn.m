function [a] = nbodyaccn(m, r)
% Solves acceleration for gravitational field.
%
%  Input arguments
%
% m:  (N x 1) Vector containing the particle masses
% r:  (N x 3) Array containing the particle positions
%
% Return argument
%
% a:  (N x 3) Array containing the computed particle accelerations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize
    N = size(r, 1);
    a = zeros(N, 3);
    M=length(m(m>0));
    G = 1;

    % Iterate to find acceleration
    for i = 1:N
        A = 0;
        for j = 1:M
            if i ~= j
                rij = [r(j, 1) - r(i, 1); 
                       r(j, 2) - r(i, 2); 
                       r(j, 3) - r(i, 3)];
                rij_mag = norm(rij);
                a_next = G * m(j) * rij / rij_mag^3;
                A = A + a_next;
            end
        end
        a(i, :) = A;
    end
end