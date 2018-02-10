% CITAMI_CALCACCELERATIONOFFSET calculates the true acceleration of the center of gravity as a function
%   of the measured accelerations and angular velocities/rates. 
% 
%              Copyright: C.C. de Visser, Delft University of Technology, 2011
%              email: c.c.devisser@tudelft.nl
%                          version 1.0
%
%   For the Citation II, the IMU has the following location:
%       x_o = 8.369 [m]  -> located at station 329.5 (8.369 [m]) 
%       y_o = 0 [m] 
%       z_o = 2.908 [m]  -> located 0.368 [m] above station 100 (2.540 [m]) 
%
function [Ax Ay Az] = citami_calcAccelerationOffset(Data, Xoffsets)

    dt = Data.tWD(2)-Data.tWD(1);
        
    Ax_ac = Data.ax;
    Ay_ac = Data.ay;
    Az_ac = Data.az;
    p = Data.p;
    q = Data.q;
    r = Data.r;
    % calculate the time derivatives of p, q, and r.
    p_dot = diff(p)*dt;
    q_dot = diff(q)*dt;
    r_dot = diff(r)*dt;
    p_dot = [p_dot; p_dot(end)];
    q_dot = [q_dot; q_dot(end)];
    r_dot = [r_dot; r_dot(end)];

    % calculate the true accelerations
    Ax = Ax_ac + Xoffsets(1).*(q.^2 + r.^2) - Xoffsets(2).*(p.*q - r_dot) - Xoffsets(3).*(p.*r + q_dot);
    Ay = Ay_ac + Xoffsets(2).*(r.^2 + p.^2) - Xoffsets(3).*(q.*r - p_dot) - Xoffsets(1).*(q.*p + r_dot);
    Az = Az_ac + Xoffsets(3).*(p.^2 + q.^2) - Xoffsets(1).*(r.*p - q_dot) - Xoffsets(2).*(r.*q + r_dot);
