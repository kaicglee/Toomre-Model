function initialize_sim1(tmax, level, num_stars, r0, v0, m)
%
% Initialization of the Toomre model simulation.
%
% Input arguments
%
% tmax:      (real scalar) Max time for the simulation to run
% level:     (real scalar) Level parameter for the FDA
% num_stars: (real scalar) # of stars per core
% v0:        (2 x 3 array) Initial velocities of cores
% m:         (2 x 1 array) Masses of cores
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initial setup
    % rotation
    cw = 1;
    
    % min and max radii of stars from core
    radmin = 0.5;
    radmax = 2;

    % generate stars
    [r01] = gen_stars(num_stars, radmin, radmax);

    % update position and velocity of stars
    r01 = r01 + r0(1, :);
    v01 = init_starmovement(m(1), r0(1, :), r01, cw);
    v01 = v01 + v0(1, :);

    m0 = [m; zeros(2*num_stars, 1)];
    r0f = [r0; r01];
    v0f = [v0; v01];

    [t, r, v] = newtongravity(tmax, level, r0f, v0f, m0);

    % reshape and plot
    hold on;
    x_1 = r(:, 1, :);
    x_1 = x_1(1, :);
    
    y_1 = r(:, 2, :);
    y_1 = y_1(1, :);

    plot(x_1, y_1, 'b');

    for i = 3:2 * num_stars + 2
        x_ = r(:, 1, :);
        x_ = x_(i, :);
        y_ = r(:, 2, :);
        y_ = y_(i, :);
        plot(x_, y_, 'r');
    end
    xlabel("x");
    ylabel("y");
    title("Star Trajectories")
    legend('Core','Star')

    function [r1] = gen_stars(num_stars, radmin, radmax)
    %
    % Generates stars randomly around both cores
    %
    % Input arguments
    %
    % num_stars (real scalar) # of stars for one core
    % radmin    (real scalar) Min radius from core
    % radmax    (real scalar) Max radius from core
    %
    % Return arguments
    %
    % r1        (num_stars x 3) Initial position vector for core 1
    % r2        (num_stars x 3) Initial position vector for core 2
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % random radius for stars around core 1 and core 2
        rad1 = (radmax - radmin) .* rand(num_stars, 1) + radmin;
        
        % random angles for stars around core 1 and core 2
        angmin = -pi;
        angmax = pi;
        ang1 = (angmax - angmin) .* rand(num_stars, 1) + angmin;

        % initialize coordinates for stars
        x = rad1 .* cos(ang1);
        y = rad1 .* sin(ang1);
        z = zeros(num_stars, 1);
        r1 = [x, y, z];
    end

    function v0 = init_starmovement(m_core, r0_core, r0_star, rotation)
    %
    % Initialize stars velocities
    %
    % Input arguments
    %
    % m:        (real scalar) Mass of core
    % r0_core:  (1 x 3) Initial position of core
    % r0_star:  (N x 3) Initial positions of stars
    % rotation: (-1 or 1) Direction of orbit
    %
    % Return arguments
    %
    % v0:       (N x 3) Initial velocities of stars
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        N = size(r0_star, 1);
        rotation = rotation * pi/2;

        % reshape
        r0_core = repmat(r0_core, N, 1);
        rij_mag = sqrt((r0_core(:,1) - r0_star(:,1)).^2 + ...
                       (r0_core(:,2) - r0_star(:,2)).^2 + ...
                       (r0_core(:,3) - r0_star(:,3)).^2);

        rz90 = [cos(rotation) -sin(rotation) 0; ...
                sin(rotation) cos(rotation) 0; ...
                0 0 1];

        sep = zeros(N, 3);
        diff = r0_core - r0_star;
        
        for j = 1:N
            sep(j, :) = (rz90 * diff(j, :).').';
        end

        v0 = sqrt(m_core ./ rij_mag) .* sep ./ rij_mag;
    end
end