function [OB_a,OT_a,WIND,PARA,take,DU,s1,s2,index] = import_data(index, alignment, add_wind, prefilter, save_data, read_data,align_method)

addpath('D:\001 - Studie\01 - Master Thesis\02_Data\OJF_november\New folder\_code_import_files');

[~,~,file] = xlsread('D:\001 - Studie\01 - Master Thesis\02_Data\OJF_november\New folder\optitrack_onboard_data_index.xlsx');

label       = file(1,:);

OT_name        = file{index,strcmp(label,'name')};
OB_name     = file{index,strcmp(label,'OB_name')}; 
wind        = file{index,strcmp(label,'wind')};
rigidbody   = file{index, strcmp(label,'rigidbody')};
configuration = file{index, strcmp(label,'configuration')};
type        = file{index,strcmp(label,'type')};
drone       = file{index,strcmp(label,'drone')};
mass        = file{index,strcmp(label,'mass')};
OT_path        = file{index,strcmp(label,'OT_path')};
OB_path        = file{index,strcmp(label,'OB_path')};
OT_function = file{index,strcmp(label,'OT_function')};
OB_function = file{index,strcmp(label,'OB_function')};
DU1          = file{index,strcmp(label,'DU1')};
DU2          = file{index,strcmp(label,'DU2')};

DU = DU1:DU2;
take.name = OT_name;
take.path = OT_path;
take.type = type;
take.windspeed = wind;
% ExtraAlign activates an extra alignment function that shifts individual
% parts 
ExtraAlign = 1;
%% read mat file if read_data is true
if read_data == 1
    load([OT_path '\Data ' OT_name,'.mat']);
else
%% import data from Optitrack and Onboard files
OT_import_function = str2func(OT_function);
OB_import_function = str2func(OB_function);

addpath(genpath('D:\001 - Studie\01 - Master Thesis\02_Data\OJF_november\New folder\_code_import_files'));
addpath(genpath('D:\001 - Studie\01 - Master Thesis\02_Data\OJF_november\_data_OJF_november_2017'));
if ~isnan(OT_name)
    addpath(OT_path); 

%     take.name_calib = '2017-11-09 05.04.16 PM';
    take.name_calib = '2016-11-21 09.57.50 AM bebop_rudi flight 2 with repeated rollpitchyaw';
    OT_filt         = OT_import_function(OT_name,take.name_calib,{rigidbody},prefilter);
    rmpath(OT_path);
    OT_a = OT_filt;
else
    OT_a = [];
end

if ~isnan(OB_name)
    addpath(OB_path);
    OB         = OB_import_function(OB_name);
    rmpath(OB_path);
    OB_a = OB;
else
    OB_a = [];
end

%% import wind speed related parameters
if wind == 'C'
    addpath(OT_path);
    T_wind = table2array(readtable(['windspeed ' OT_name '.txt'],'Delimiter','\t'));
    time_wind = T_wind(:,1);
    speed_wind = T_wind(:,5);
    time_wind(time_wind==0) = []; speed_wind(speed_wind==0) = [];

    take_date = str2double(OT_name(9:10));
    take_hour = str2double(OT_name(12:13));
    take_minute = str2double(OT_name(15:16));
    take_sec = str2double(OT_name(18:19));
    take_AMPM = OT_name(21:22);
    if take_AMPM == 'PM'
        take_hour = take_hour + 12;
    end
    % relative time to 2017-11-08 13.00.00 PM
    take_time = (take_date - 8)*3600*24 + (take_hour - 13)*3600 + take_minute*60 + take_sec; 

    WIND.speed = speed_wind;
    WIND.TIME = time_wind-take_time; %align with OT
    rmpath(OT_path);
else
    WIND.TIME = OT_filt.TIME;
    WIND.speed = wind*ones(size(OT_filt.TIME));
end
%% resample and finddelay between OB and OT
if alignment == 1
    OT = struct;
    allnames = fieldnames(OT_filt);
    OB_discard = find(OB.TIME > OT_filt.TIME(end));
    time_OB = OB.TIME;
    time_OB(OB_discard) = [];
    for k = 1:length(fieldnames(OT_filt))
        ct_orig = OT_filt.(allnames{k});
        ct_orig = interpNan(ct_orig);
        OT.(allnames{k}) = interp1(OT_filt.TIME,ct_orig,time_OB,'spline');
    end

    switch configuration
        case 'nominal'
%             method = 'theta';
            method = align_method;
        otherwise
            method = 'phi_ot';
    end
    [OT_a,OB_a,Delay,s1,s2]                     = align_signals(OT,OB,method);
    
    if (ExtraAlign == 1) 
    [OT_shift,OB_shift,s1,s2,shifting_vector]   = align_signal2(OT_a,OB_a,method);
    OT_a = OT_shift;
    OB_a = OB_shift; 
    s1 = OT_a.PHI; s2 = OB_a.phi_ot*57.3;
    end
end

%% resemble WIND and align with OT_a
if add_wind == 1 && sum(~isnan(OT_name)) 
    
    if Delay >=0
        %align with OB (otherwise align with OT, WIND.TIME unchange)
         WIND.TIME = WIND.TIME - Delay/512;
    end
    
%     speed_wind_resemple = interp1(WIND.TIME,WIND.speed,OT_a.TIME,'spline');
    speed_wind_resemple = LinearInterpWithClipExtrap(WIND.TIME,WIND.speed,OT_a.TIME);
    OT_a.VY_air = OT_a.VY - speed_wind_resemple;
    OT_a.VX_air = OT_a.VX;
    OT_a.VZ_air = OT_a.VZ;
    
    VwO = zeros(length(OT_a.TIME),3);
    VwB = zeros(length(OT_a.TIME),3);
    for i = 1:length(OT_a.TIME)
        Rot_G2O = Rot_z(-OT_a.att_G2O(i,1)*pi/180)*Rot_x(-OT_a.att_G2O(i,2)*pi/180)*Rot_y(-OT_a.att_G2O(i,3)*pi/180);
        Vw = [speed_wind_resemple(i) 0 0]';
        VwO(i,:) = (Rot_G2O*Vw)';
    end
    VwB(:,1) = VwO(:,3);
    VwB(:,2) = -VwO(:,1);
    VwB(:,3) = -VwO(:,2);
        
    OT_a.U_air = OT_a.U + VwB(:,1);
    OT_a.V_air = OT_a.V + VwB(:,2);
    OT_a.W_air = OT_a.W + VwB(:,3);
else
    OT_a.VX_air = OT_a.VX;
    OT_a.VY_air = OT_a.VY;
    OT_a.VZ_air = OT_a.VZ;
    OT_a.U_air = OT_a.U;
    OT_a.V_air = OT_a.V;
    OT_a.W_air = OT_a.W;
end

end
%% import drone geometry parameters
PARA.mass = mass;
addpath('D:\001 - Studie\01 - Master Thesis\02_Data\OJF_november\New folder\_data_drone_parameters');
switch type
    case 'Bebop'
        load('TUBB20_bare_parameters.mat');
    case 'Bebop2'
        load('Bebop2_leon_parameters.mat');
%         load('TUBB20_parameters.mat');
    otherwise
end
rmpath('E:\Data\_data_drone_parameters');

PARA.Iv = parameters.Iv;
PARA.b  = parameters.b;
PARA.l  = parameters.l;
PARA.Ip = parameters.Ip;
PARA.R  = parameters.R;

%% save data into a certain path
if save_data == true && read_data == false
    save([OT_path '\Data ' OT_name,'.mat'],'OT_a','OB_a','WIND','PARA');
end
end


function [OT_a,OB_a,Delay,s1,s2] = align_signals(OT,OB,method)

switch method
    case 'phi_ot'
        [s1,s2,Delay] = alignsignals(OB.phi_ot*57.3,OT.PHI);
    case 'theta_ot'
        [s1,s2,Delay] = alignsignals(OB.theta_ot*57.3,OT.THETA);        
    case 'phi'
        [s1,s2,Delay] = alignsignals(OB.phi*57.3,OT.PHI);
    case 'theta'
        [s1,s2,Delay] = alignsignals(OB.theta*57.3,OT.THETA);          
    case 'Z'
        [s1,s2,Delay] = alignsignals(OB.PosNED(:,3),-OT.posCO_G(:,2)-(-OT.posCO_G(1,2)-OB.PosNED(1,3))); %using height
    case 'Z_top' %Use the upper 20% height to align signals
        hthr = min(OB.PosNED(:,3))*0.8;
        h1 = find(OB.PosNED(:,3)<=hthr); h2 = find(-OT.posCO_G(:,2)-(-OT.posCO_G(1,2)-OB.PosNED(1,3))<=hthr);
        hh1 = -0*ones(size(OB.PosNED(:,3))); hh2 = -0*ones(size(-OT.posCO_G(:,2)));
        hh1(h1) = OB.PosNED(h1,3); hh2(h2) = -OT.posCO_G(h2,2)-(-OT.posCO_G(1,2)-OB.PosNED(1,3));
        [s1,s2,Delay] = alignsignals(-hh1,-hh2);
    case 'PQR'        
        pqr_filt_ot = [butterworth(OT.P,4,5/512);butterworth(OT.Q,4,5/512);butterworth(OT.R,4,5/512)];
        pqr_filt_ob = [butterworth(OB.p,4,5/512);butterworth(OB.q,4,5/512);butterworth(OB.r,4,5/512)];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot);
    case 'P'        
        pqr_filt_ot = [butterworth(OT.P,4,5/512);];
        pqr_filt_ob = [butterworth(OB.p,4,5/512);];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot);
    case 'Q'        
        pqr_filt_ot = [butterworth(OT.Q,4,5/512);];
        pqr_filt_ob = [butterworth(OB.q,4,5/512);];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot); 
    case 'R'        
        pqr_filt_ot = [butterworth(OT.R,4,5/512);];
        pqr_filt_ob = [butterworth(OB.r,4,5/512);];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot);         
    case 'VZPQR'
        pqr_filt_ot = [butterworth(OT.P,4,5/512);butterworth(OT.Q,4,5/512);butterworth(OT.R,4,5/512)];
        pqr_filt_ob = [butterworth(OB.P,4,5/512);butterworth(OB.Q,4,5/512);butterworth(OB.R,4,5/512)];
        [s1,s2,Delay] = alignsignals([pqr_filt_ob;OB.VelNED(:,3)],...
                                     [pqr_filt_ot*57.3;OT.vel_E(:,3)]);
end
% %%
%%
L_OB = length(OB.TIME);
L_OT  = length(OT.TIME);

fields_OT = fieldnames(OT);
fields_OB = fieldnames(OB);

OT_a = struct; OB_a = struct;
if Delay >= 0 
    % data aligned to OB
    L_SYN = min(L_OT-Delay,L_OB);
    for i = 1:length(fields_OT)
        OT_a.(fields_OT{i}) = OT.(fields_OT{i})(Delay+1:Delay+L_SYN,:);
    end
    for i = 1:length(fields_OB)
        OB_a.(fields_OB{i}) = OB.(fields_OB{i})(1:L_SYN,:);
    end
    OT_a.TIME = OB_a.TIME;
else   
    % data aligned to OT
    L_SYN = min(L_OT,L_OB+Delay);
    for i = 1:length(fields_OT)
        OT_a.(fields_OT{i}) = OT.(fields_OT{i})(1:L_SYN,:);
    end
    for i = 1:length(fields_OB)
        OB_a.(fields_OB{i}) = OB.(fields_OB{i})(1-Delay:L_SYN-Delay,:);
    end
    OB_a.TIME = OT_a.TIME;
end
%%
% Check align result
figure('position',[0 0,700,200])
% plot(OB_a.PosNED(:,3)); hold on
% plot(-OT_a.posCO_G(:,2)-(-OT_a.posCO_G(1,2)-OB_a.PosNED(1,3))); ylabel(method);title('Check alignment');
plot(s1);hold on;
plot(s2);
end

function vq = LinearInterpWithClipExtrap(x,v,xq)

    vq = interp1(x,v,xq);

    [XMax, idxVMax] = max(x);
    [XMin, idxVMin] = min(x);

    idxMax = xq > XMax;
    idxMin = xq < XMin;

    vq(idxMax) = v(idxVMax);
    vq(idxMin) = v(idxVMin);

end