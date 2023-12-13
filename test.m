% script to run test for newtongravity for two particles in circular orbit

% mass
m1 = 1;
m2 = 0.5;
m_list = [m1; m2];
m = sum(m_list);

% radius
r = 0.75;
r1 = m2 / m * r;
r2 = m1 / m * r;
r0 = [[r1, 0, 0]; [-r2, 0, 0]];

% velocity
v1 = sqrt(m2 * r1) / r;
v2 = sqrt(m1 * r2) / r;
v0 = [[0, v1, 0]; [0, -v2, 0]];

% time
tmax = 40;
level = 8;
nt = 2^level + 1;

% Call function
[t,r,v]=newtongravity(tmax, 8, r0, v0, m_list);

% Reshape first particle
x1 = r(:, 1,  1:end);
x1 = x1(1, 1:end);

y1 = r(:, 2,  1:end);
y1 = y1(1, 1:end);

% Reshape second particle
x2 = r(:, 1,  1:end);
x2 = x2(2, 1:end);

y2 = r(:, 2,  1:end);
y2 = y2(2, 1:end);

% Plot x and y coordinates of two particles
clf;
hold on;
xlabel("x");
ylabel("y");
plot(x1, y1);
plot(x2, y2);
legend('Particle 1 (m = 1.0)', 'Particle 2 (m = 0.5)');
