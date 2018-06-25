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

OB.p = data(:,31);
OB.q = data(:,32);
OB.r = data(:,33);

OB.phi = data(:,34);
OB.theta = data(:,35);
OB.psi = data(:,36);

OB.phi_ot = data(:,40);
OB.theta_ot = data(:,41);
OB.psi_ot = data(:,42);
OB.r_ot = data(:,43);
OB.p_des_pa = data(:,44);
OB.q_des_pa = data(:,45);
OB.r_des_pa = data(:,46);

OB.h1 = data(:,47);
OB.h2 = data(:,48);
OB.h3 = data(:,49);

OB.ndix = data(:,50);
OB.ndiy = data(:,51);
OB.ndiz = data(:,52);

OB.ax_des = data(:,53);
OB.ay_des = data(:,54);
OB.az_des = data(:,55);

OB.ax_des_filter = data(:,56);
OB.ay_des_filter = data(:,57);
OB.az_des_filter = data(:,58);

OB.p_dot_des = data(:,59);
OB.q_dot_des = data(:,60);
OB.r_dot_des = data(:,61);

OB.p_des = data(:,62);
OB.q_des = data(:,63);
OB.r_des = data(:,64);

OB.z_ref = data(:,65);
OB.z 	 = data(:,66);

OB.du1 = data(:,67);
OB.du2 = data(:,68);
OB.du3 = data(:,69);
OB.du4 = data(:,70);

OB.Vx = data(:,71);
OB.Vy = data(:,72);
OB.Vz = data(:,73);
OB.Vx_des = data(:,74);
OB.Vy_des = data(:,75);
OB.Vz_des = data(:,76);

OB.x_ref = data(:,77);
OB.y_ref = data(:,78);
OB.x = data(:,79);
OB.y = data(:,80);

OB.forward_switch= data(:,83);
end
