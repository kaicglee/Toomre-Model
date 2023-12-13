% script to run the convergence test for newtongravity

% time
tmax = 100;

% mass
m1 = 1;
m2 = 0.5;
m_list = [m1; m2];
m = sum(m_list);

% radius
r = 4;
r1 = m2 / m * r;
r2 = m1 / m * r;
r0 = [[r1, 0, 0];[-r2, 0, 0]];

% velocity
v1 = sqrt(m2 * r1) / r;
v2 = sqrt(m1 * r2) / r;
v0 = [[0, v1, 0];[0, -v2, 0]];

% Call function for each level
[t6 r6 v6]   =  newtongravity(tmax, 6, r0, v0, m_list);
[t7 r7 v7]   =  newtongravity(tmax, 7, r0, v0, m_list);
[t8 r8 v8]   =  newtongravity(tmax, 8, r0, v0, m_list);
[t9 r9 v9]   =  newtongravity(tmax, 9, r0, v0, m_list);

% Reshape to get only x-coordinates
r6 = r6(:, 1,  1:end);
r7 = r7(:, 1,  1:end);
r8 = r8(:, 1,  1:end);
r9 = r9(:, 1,  1:end);

r6 = r6(1, 1:end);
r7 = r7(1, 1:end);
r8 = r8(1, 1:end);
r9 = r9(1, 1:end);

% Plot convergence test for three levels
clf; 
figure(1);
hold on; 
plot(t6, r6, 'r-.o');
plot(t7, r7, 'g-.+'); 
plot(t8, r8, 'b-.*');
xlabel("Time");
ylabel("x");
legend('level 6', 'level 7', 'level 8')

% Reshape to get same sized vectors
r7 = r7(1:2:end);
r8 = r8(1:4:end);
r9 = r9(1:8:end);

% Take differences between levels
dr67 = r6 - r7;
dr78 = r7 - r8;
dr89 = r8 - r9;

% Plot convergence test for differences
figure(2)
hold on; 
plot(t6, dr67, 'r-.o'); 
plot(t6, dr78, 'g-.+');
plot(t6, dr89, 'b-.*'); 
xlabel("Time");
legend('level 6-7', 'level 7-8', 'level 8-9')

% Scale differences
dr78 = 4 * dr78;
dr89 = 16 * dr89;

% Plot conveergence test for scaled differences
figure(3)
hold on; 
plot(t6, dr67, 'r-.o'); 
plot(t6, dr78, 'g-.+');
plot(t6, dr89, 'b-.*'); 
xlabel("Time");
legend('level 6-7', '4*level 7-8', '16*level 8-9')