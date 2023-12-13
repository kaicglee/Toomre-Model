function [t r v] = newtongravity(tmax, level, r0, v0, m)
%  Solves Toomre Model using second order finite difference approximation.
%
%  Input arguments
%
%      tmax:  (real scalar) Final solution time.
%      level: (integer scalar) Discretization level.
%      r0:    (real vector) Initial positions of the particles.
%      v0:    (real vector) Initial velocities of the particles.
%
%  Return values
%
%      t:     (real vector) Vector of length nt = 2^level + 1 containing
%             discrete times (time mesh).
%      r:     (real vector) Vector of length (N, 3, nt) containing computed
%             positions at discrete times t(n).
%      v:     (real vector) Vector of length (N, 3, nt) containing computed
%             velocities at discrete times t(n).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % trace controls "tracing" output.  Set 0 to disable, non-0 to enable.
   trace = 1;
   % tracefreq controls frequency of tracing output in main time step loop.
   tracefreq = 100;

   if trace
      fprintf('In newtongravity: Argument dump follows\n');
      tmax, level, r0, v0, m
   end

   % Define number of time steps and create t, r and v arrays of
   % appropriate size.
   % N:   number of particles
   % nt:  total number of time steps
   nt = 2^level + 1;
   t = linspace(0.0, tmax, nt);
   N = size(r0, 1);
   r = zeros(N, 3, nt);
   v = zeros(N, 3, nt);

   % Determine discrete time step from t array.
   deltat = t(2) - t(1);

   % Initialize first values of the particles accelerations.
   a0 = nbodyaccn(m, r0);

   % Initialize first values of the particles positions.
  
   r(:, :, 1) = r0;
   r(:, :, 2) = r0 + deltat * v0 + 0.5 * deltat^2 * a0;

   if trace
      fprintf('deltat=%g r(:, :, 1)=%g r(:, :, 2)=%g\n',...
              deltat, r(:, :, 1), r(:, :, 2));
   end

   % Initialize first values of the particles velocity.
   v(:, :, 1) = v0;

   % Evolve the Toomre model to the final time using the discrete equations
   % of motion.  Also compute an estimate of the velocity at
   % each time step.

   for n = 2 : nt-1
      % This generates tracing output every 'tracefreq' steps.
      if rem(n, tracefreq) == 0
         fprintf('nbody: Step %d of %d\n', n, nt);
      end

      a_n = nbodyaccn(m, r(:, :, n));
      
      r(:, :, n+1) = 2 * r(:, :, n) - r(:, :, n-1) + deltat^2 * a_n;
      v(:, :, n) = (r(:, :, n+1) - r(:, :, n-1)) / (2 * deltat);
   end

   % Use linear extrapolation to determine the value of r at the
   % final time step.

   v(:, :, nt) = 2 * v(:, :, nt-1) - v(:, :, nt-2);

end