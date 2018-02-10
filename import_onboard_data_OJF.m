function OB = import_onboard_data_OJF(name)

OB = struct;

data = csvread([name '.csv'],1,0);
% 
S_acc = 4096;
S_om = 1879.3; 


OB.TIME = (1:size(data,1))'/512;

OB.p = data(:,2)/S_om;
OB.q = data(:,3)/S_om;
OB.r = data(:,4)/S_om;
OB.P = data(:,2)/S_om;
OB.Q = data(:,3)/S_om;
OB.R = data(:,4)/S_om;
OB.ax = data(:,5)/S_acc;
OB.ay = data(:,6)/S_acc;
OB.az = data(:,7)/S_acc;
OB.AX = data(:,5)/S_acc;
OB.AY = data(:,6)/S_acc;
OB.AZ = data(:,7)/S_acc;

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

OB.w1obs_indi = data(:,27);
OB.w2obs_indi = data(:,28);
OB.w3obs_indi = data(:,29);
OB.w4obs_indi = data(:,30);

end