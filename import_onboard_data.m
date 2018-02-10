function OB = import_onboard_data(name)

OB = struct;

data = csvread(['onboard_log ' name '.csv'],0,0);
% IMU = struct;
% OTOB = struct;
% BAT = struct;
% COM = struct;
% ACT = struct;
g = 9.8124;
scale_pqr = 0.0139882;   %deg/s
scale_acc = 1/1024/g; %1/1024

S_acc = 4096; % 2048 according to datasheet
S_om = 1879.3; % 16.4 according to datasheet

OB.COUNTER     = data(:,1);
OB.TIME    = data(:,2);
OB.P       = data(:,4) * scale_pqr;
OB.Q       = data(:,5) * scale_pqr;
OB.R       = data(:,6) * scale_pqr;
OB.AX      = data(:,7) * scale_acc;
OB.AY      = data(:,8) * scale_acc;
OB.AZ      = data(:,9) * scale_acc;
% OB.p = data(:,10) / S_om;
% OB.q = data(:,11) / S_om;
% OB.r = data(:,12) / S_om;
% OB.ax = data(:,13) / S_acc;
% OB.ay = data(:,14) / S_acc;
% OB.az = data(:,15) / S_acc;
OB.MAGX    = data(:,16);
OB.MAGY    = data(:,17);
OB.MAGZ    = data(:,18);
OB.phi     = data(:,19)*57.3;
OB.theta   = data(:,20)*57.3;
OB.psi     = data(:,21)*57.3;

OB.PosNED = data(:,22:24);
OB.VelNED = data(:,25:27);

% OB.COM.THRUST = data(:,28);
% OB.COM.PHI = data(:,29);
% OB.COM.THETA = data(:,30);
% OB.COM.PSI = data(:,31);
% OB.BAT.VOLTAGE = data(:,32);
% OB.BAT.CURRENT = data(:,33);
OB.w1ref = data(:,34);
OB.w2ref = data(:,35);
OB.w3ref = data(:,36);
OB.w4ref = data(:,37);

OB.RPM1 = data(:,38);
OB.RPM2 = data(:,39);
OB.RPM3 = data(:,40);
OB.RPM4 = data(:,41);
end