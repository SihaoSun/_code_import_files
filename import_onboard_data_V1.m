function OB = import_onboard_data_V1(name)

% first version onboard data importor with no damage controller related
% variables logged.

OB = struct;

data = csvread([name '.csv'],1,0);
% IMU = struct;
% OTOB = struct;
% BAT = struct;
% COM = struct;
% ACT = struct;
% g = 9.8124;
% scale_pqr = 0.0139882;   %deg/s
% scale_acc = 1/1024/g; %1/1024
% 
% S_acc = 4096; % 2048 according to datasheet
% S_om = 1879.3; % 16.4 according to datasheet


OB.TIME = (1:size(data,1))'/512;

OB.thrcmd = data(:,11);

OB.w1obs = data(:,19);
OB.w2obs = data(:,20);
OB.w3obs = data(:,21);
OB.w4obs = data(:,22);

OB.w1ref = data(:,23);
OB.w2ref = data(:,24);
OB.w3ref = data(:,25);
OB.w4ref = data(:,26);

OB.w1obs_indi = data(:,27);
OB.w2obs_indi = data(:,28);
OB.w3obs_indi = data(:,29);
OB.w4obs_indi = data(:,30);

OB.p = data(:,31);
OB.q = data(:,32);
OB.r = data(:,33);

OB.phi = data(:,34);
OB.theta = data(:,35);
OB.psi = data(:,36);

OB.ax = data(:,37)/1000;
OB.ay = data(:,38)/1000;
OB.az = data(:,39)/1000;

OB.phi_ot = data(:,40);
OB.theta_ot = data(:,41);
OB.psi_ot = data(:,42);
OB.r_ot = data(:,43);

OB.p_des = data(:,44);
OB.q_des = data(:,45);
OB.r_des = data(:,46);

% OB.h1 = data(:,47);
% OB.h2 = data(:,48);
% OB.ndix = data(:,49);
% OB.ndiy = data(:,50);
% OB.ndiz = data(:,51);
% OB.acc_des_x = data(:,52);
% OB.acc_des_y = data(:,53);
% OB.acc_des_z = data(:,54);
% OB.acc_des_x_filter = data(:,55);
% OB.acc_des_y_filter = data(:,56);
% OB.acc_des_z_filter = data(:,57);
% OB.p_des_dot = data(:,58);
% OB.q_des_dot = data(:,59);
% OB.p_des_filter = data(:,60);
% OB.q_des_filter = data(:,61);

OB.P = OB.p;
OB.Q = OB.q;
OB.R = OB.r;
OB.AX = OB.ax;
OB.AY = OB.ay;
OB.AZ = OB.az;
end