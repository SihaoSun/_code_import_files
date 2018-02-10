function [OB] = import_onboard_data_rudi(take,string)
if nargin == 0
    take = select_take(20);
    do_filtering = true;
elseif nargin == 1 %default
    do_filtering = false;
elseif nargin == 2
    if strcmp(string,'filter_on')
        do_filtering = true;
    else
        do_filtering = false;
    end
end
% take = select_take(12)
% do_filtering = true;
%% Script Parameters

IMU = struct;
INT32_RATE_FRAC = 2^-12; % converting int32 to float for gyro
INT32_ACCEL_FRAC = 2^-10; % converting int32 to float for accelerometer
g = 9.8124;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import on-board data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load on-board data
logdata = csvread(['onboard_log ' take.name '.csv'],1,0);

% Retrieving labels of the logged data from first line
openedFile = fopen(['onboard_log ' take.name '.csv'],'r');
line = fgetl(openedFile);
logdata_label = regexp(line, ',', 'split');
fclose(openedFile);

% Parameters and calibration constants
g = 9.8124;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter corrupt lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) remove corrupt lines at the end starting with the first line full of zeros
i = 1;
ilast = 0;
while ilast == 0
    if logdata(i,:)==0
        ilast = i-1;
    end
    i = i+1;
    if i>length(logdata) % no lines full of zeros
        ilast = -1;
    end
end

if ilast==-1
    index = [1:length(logdata)]';
else
    index = [1:ilast]';
end

IMU.time = (logdata(index,strcmp(logdata_label,'time_us'))-logdata(1,strcmp(logdata_label,'time_us')))*1e-6; % transformed to seconds 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) remove lines with incorrect time stamp

dt = [diff(IMU.time); 1]; % positive value added at the end (this would be the case if another correct timestamp was present)
lines2check=index(dt<0);
for i = 1:length(lines2check)
    if lines2check(i)>1 && lines2check(i)<length(index) && dt(lines2check(i)-1)<dt(lines2check(i)+1) 
    % dt(i-1)<dt(i+1)  -->  time is lower than it should be - the sample to be removed is i+1
        lines2check(i) = lines2check(i)+1;
    else    
    % dt(i-1)>dt(i+1)  -->  time is higher than it should be - no change
        lines2check(i) = lines2check(i);
    end
end

lines2remove = index(dt>1)+1; % lines where timestamp increased by at least one second
    
j = 1;
i = 1;
index2 = [];
while i<=length(index)
    if j<=length(lines2check) && lines2check(j)==index(i)
        % skip lines with negative dt
        j = j+1;
        i = i+1;
    elseif ~isempty(lines2remove) && index(i)>=lines2remove(1)
        % skip lines after too positive dt
        i = i+1;
    else
        index2 = [index2; index(i)];
        i = i+1;
    end
end

index = index2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve data after filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSample     = (logdata(index,1));

IMU.time    = (logdata(index,strcmp(logdata_label,'time_us'))-logdata(1,strcmp(logdata_label,'time_us')))*1e-6; % transformed to seconds 

ACC.x_raw   = (logdata(index,strcmp(logdata_label,'accel_unscaled_x'))); % [g]
ACC.y_raw   = (logdata(index,strcmp(logdata_label,'accel_unscaled_y'))); % [g]
ACC.z_raw   = (logdata(index,strcmp(logdata_label,'accel_unscaled_z'))); % [g]

GYR.p_raw   = (logdata(index,strcmp(logdata_label,'gyro_unscaled_p'))); % [rad/s]
GYR.q_raw   = (logdata(index,strcmp(logdata_label,'gyro_unscaled_q'))); % [rad/s]
GYR.r_raw   = (logdata(index,strcmp(logdata_label,'gyro_unscaled_r'))); % [rad/s]

ACC.x   = (logdata(index,strcmp(logdata_label,'accel_x')))*INT32_ACCEL_FRAC; % [m/s^2]
ACC.y   = (logdata(index,strcmp(logdata_label,'accel_y')))*INT32_ACCEL_FRAC; % [m/s^2]
ACC.z   = (logdata(index,strcmp(logdata_label,'accel_z')))*INT32_ACCEL_FRAC; % [m/s^2]

GYR.p   = (logdata(index,strcmp(logdata_label,'gyro_p')))*INT32_RATE_FRAC; % [rad/s]
GYR.q   = (logdata(index,strcmp(logdata_label,'gyro_q')))*INT32_RATE_FRAC; % [rad/s]
GYR.r   = (logdata(index,strcmp(logdata_label,'gyro_r')))*INT32_RATE_FRAC; % [rad/s]

ACC.N_state   = (logdata(index,strcmp(logdata_label,'accel_x_state'))); %floats [m/s^2]
ACC.E_state   = (logdata(index,strcmp(logdata_label,'accel_y_state'))); %floats [m/s^2]
ACC.D_state   = (logdata(index,strcmp(logdata_label,'accel_z_state'))); %floats [m/s^2]

MAG.x   = (logdata(index,strcmp(logdata_label,'mag_unscaled_x')));
MAG.y   = (logdata(index,strcmp(logdata_label,'mag_unscaled_y')));
MAG.z   = (logdata(index,strcmp(logdata_label,'mag_unscaled_z')));

B2IMU.phi   = (logdata(index,strcmp(logdata_label,'body2imu_x'))); % [?]
B2IMU.theta = (logdata(index,strcmp(logdata_label,'body2imu_y'))); % [?]
B2IMU.psi   = (logdata(index,strcmp(logdata_label,'body2imu_z'))); % [?]

CMD.T       = (logdata(index,strcmp(logdata_label,'COMMAND_THRUST')));  % [?]
CMD.roll    = (logdata(index,strcmp(logdata_label,'COMMAND_ROLL'))); % [?]
CMD.pitch   = (logdata(index,strcmp(logdata_label,'COMMAND_PITCH'))); % [?]
CMD.yaw     = (logdata(index,strcmp(logdata_label,'COMMAND_YAW'))); % [?]

ATT.qw      = (logdata(index,strcmp(logdata_label,'qi'))); % [-]
ATT.qx      = (logdata(index,strcmp(logdata_label,'qx'))); % [-]
ATT.qy      = (logdata(index,strcmp(logdata_label,'qy'))); % [-]
ATT.qz      = (logdata(index,strcmp(logdata_label,'qz'))); % [-]

ROTOR.w1    = (logdata(index,strcmp(logdata_label,'w1'))); % [rpm]
ROTOR.w2    = (logdata(index,strcmp(logdata_label,'w2'))); % [rpm]
ROTOR.w3    = (logdata(index,strcmp(logdata_label,'w3'))); % [rpm]
ROTOR.w4    = (logdata(index,strcmp(logdata_label,'w4'))); % [rpm]

POS.x       = (logdata(index,strcmp(logdata_label,'state_posX'))); %[m];
POS.y       = (logdata(index,strcmp(logdata_label,'state_posY'))); %[m];
POS.z       = (logdata(index,strcmp(logdata_label,'state_posZ'))); %[m];

SPEED.x     = (logdata(index,strcmp(logdata_label,'state_speed_N'))); %[m];
SPEED.y     = (logdata(index,strcmp(logdata_label,'state_speed_E'))); %[m];
SPEED.z     = (logdata(index,strcmp(logdata_label,'state_speed_D'))); %[m];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adjust for paparazzi orientation, is to be: xyz = forward, right, down
axg_raw         = ACC.x_raw;
ayg_raw         = ACC.y_raw;
azg_raw         = ACC.z_raw;
omx_raw         = GYR.p_raw;
omy_raw         = GYR.q_raw;
omz_raw         = GYR.r_raw;

% Mean sampling frequency
fs = round(mean(1./diff(IMU.time)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Converting units from unscaled data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion bit --> physical units
%S_acc = sqrt(mean(ax)^2+mean(ay)^2+mean(az)^2) - only gravity is measured in still
S_acc = 4096; % 4096 (datasheet)
S_om = 32.8*180/pi; %1879.3 or 4096/2.17953 or 32.8*180/pi (datasheet, sensitivy level 2)

axg = axg_raw/S_acc; % in g
ayg = ayg_raw/S_acc;
azg = azg_raw/S_acc;
omx = omx_raw/S_om; % in rad/s
omy = omy_raw/S_om;
omz = omz_raw/S_om;

%% Filter accelerations if desired
if do_filtering
    [b,a] = butter(4,2*40/512);
    axg = filtfilt(b,a,axg);
    ayg = filtfilt(b,a,ayg);
    azg = filtfilt(b,a,azg);
%     omx = filtfilt(b,a,omx);
%     omy = filtfilt(b,a,omy);
%     omz = filtfilt(b,a,omz);
end

r2d = 180/pi;

OB.COUNTER     = nSample;
OB.TIME    = IMU.time;
OB.IMU.TIME    = IMU.time;
OB.IMU.P       = GYR.p * r2d;
OB.IMU.Q       = GYR.q * r2d;
OB.IMU.R       = GYR.r * r2d;
OB.IMU.AX      = axg;
OB.IMU.AY      = ayg;
OB.IMU.AZ      = azg;

% OB.IMU.PHI     = data(:,19)*57.3;
% OB.IMU.THETA   = data(:,20)*57.3;
% OB.IMU.PSI     = data(:,21)*57.3;

OB.OTOB.PosNED = [POS.x POS.y POS.z];
OB.OTOB.VelNED = [SPEED.x SPEED.y SPEED.z];

% OB.COM.THRUST = data(:,28);
% OB.COM.PHI = data(:,29);
% OB.COM.THETA = data(:,30);
% OB.COM.PSI = data(:,31);
% OB.BAT.VOLTAGE = data(:,32);
% OB.BAT.CURRENT = data(:,33);
% OB.COM.RPM1 = data(:,34);
% OB.COM.RPM2 = data(:,35);
% OB.COM.RPM3 = data(:,36);
% OB.COM.RPM4 = data(:,37);

OB.ACT.RPM1 = ROTOR.w1;
OB.ACT.RPM2 = ROTOR.w2;
OB.ACT.RPM3 = ROTOR.w3;
OB.ACT.RPM4 = ROTOR.w4;

return
end