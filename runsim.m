% script to run initialize_sim

tmax = 40;
level = 8;
num_stars = 3000;

m1 = 10;
m2 = 4;
m_list = [m1; m2];

r0 = [[2, 2, 0]; [10, 2, 0]];
v0 = [[0, 0, 0]; [-0.5, 0, 0]];

initialize_sim(tmax, level, num_stars, r0, v0, m_list);