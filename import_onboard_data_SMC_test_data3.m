function OB = import_onboard_data_SMC_test_data3(name)

OB = struct;

data = csvread([name '.csv'],1,0);
% 
S_acc = 4096;
S_om = 1879.3; 


OB.TIME = (1:size(data,1))'/512;

OB.p = data(:,2)/S_om;
OB.q = data(:,3)/S_om;
OB.r = data(:,4)/S_om;
OB.ax = data(:,5)/S_acc;
OB.ay = data(:,6)/S_acc;
OB.az = data(:,7)/S_acc;
OB.CMD_THR = data(:,11);
OB.CMD_ROLL = data(:,12);
OB.CMD_PITCH = data(:,13);
OB.CMD_YAW = data(:,14);

OB.qi = data(:,15);
OB.qx = data(:,16);
OB.qy = data(:,17);
OB.qz = data(:,18);

Q = data(:,15:18);
[psi, theta, phi] = quat2angle(Q);

OB.psi = psi;
OB.theta = theta;
OB.phi = phi;

OB.w1obs = data(:,19);
OB.w2obs = data(:,20);
OB.w3obs = data(:,21);
OB.w4obs = data(:,22);

OB.w1ref = data(:,23);
OB.w2ref = data(:,24);
OB.w3ref = data(:,25);
OB.w4ref = data(:,26);

OB.act_obs1 = data(:,27);
OB.act_obs2 = data(:,28);
OB.act_obs3 = data(:,29);
OB.act_obs4 = data(:,30);

OB.phi = data(:,34);
OB.theta = data(:,35);
OB.psi = data(:,36);

OB.phi_ot = data(:,40);
OB.theta_ot = data(:,41);
OB.psi_ot = data(:,42);
OB.h1 = data(:,47);
OB.h2 = data(:,48);
OB.ndix = data(:,49);
OB.ndiy = data(:,50);
OB.ndiz = data(:,51);

OB.p_dot_des = data(:,58);
OB.q_dot_des = data(:,59);
OB.r_dot_des = data(:,60);

OB.p_des = data(:,61);
OB.q_des = data(:,62);
OB.r_des = data(:,63);

OB.nu01  = data(:,64);
OB.nu02  = data(:,65);
OB.nu03  = data(:,66);
OB.nu04  = data(:,67);

OB.nueq1 = data(:,68);
OB.nueq2 = data(:,69);
OB.nueq3 = data(:,70);
OB.nueq4 = data(:,71);

OB.PSI01 = data(:,72);
OB.PSI02 = data(:,73);
OB.PSI03 = data(:,74);
OB.PSI04 = data(:,75);

OB.zdot1 = data(:,76);
OB.zdot2 = data(:,77);
OB.zdot3 = data(:,78);
OB.zdot4 = data(:,79);

OB.sigma1 = data(:,80);
OB.sigma2 = data(:,81);
OB.sigma3 = data(:,82);
OB.sigma4 = data(:,83);

OB.s1 = data(:,84);
OB.s2 = data(:,85);
OB.s3 = data(:,86);
OB.s4 = data(:,87);

% % OB.step = data(:,88);
OB.nx_des_step = data(:,88);
OB.ny_des_step = data(:,89);

OB.e1 = data(:,90);
OB.e2 = data(:,91);
OB.e3 = data(:,92);
OB.e4 = data(:,93);

OB.z_ref = data(:,94);
OB.z = data(:,95);

OB.du1 = data(:,96);
OB.du2 = data(:,97);
OB.du3 = data(:,98);
OB.du4 = data(:,99);

if str2double(name) >= 11
OB.Vx = data(:,100);
OB.Vy = data(:,101);
OB.Vz = data(:,102);
OB.Vx_des = data(:,103);
OB.Vy_des = data(:,104);
OB.Vz_des = data(:,105);
end
end