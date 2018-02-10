%% VERIFICATION
if (VERIFICATION == 1)
    % Take the opposing state to cross check if the shift went okay
switch method
    case 'phi_ot'
        veri_1 = OB_a.theta_ot*57.3;
        veri_2 = OT_a.THETA;
        maxlag_firstloop = 300;
        maxlag_otherloops = 100;
    case 'theta_ot'

        veri_1 = OB_a.phi_ot*57.3;
        veri_2 = OT_a.PHI;
        maxlag_firstloop = 300;
        maxlag_otherloops = 100;
    case 'phi'
        veri_1 = OB_a.theta*57.3;
        veri_2 = OT_a.THETA;
        maxlag_firstloop = 3;
        maxlag_otherloops = 1;
    case 'theta'
        veri_1 = OB_a.phi*57.3;
        veri_2 = OT_a.PHI;
        maxlag_firstloop = 3;
        maxlag_otherloops = 1;
end 

figure; hold on; plot(veri_1); plot(veri_2); 

veri_11 = veri_1(crop_vector(1,1):crop_vector(2,1));
veri_22 = veri_2(crop_vector(1,1):crop_vector(2,1));
plot(veri_11); plot(veri_22); hold off;

% Perform first shift outside loop
if (shifting_vector(2,1) < 0)
    % OB should shift to the right, thus inserting zeros at the start
    % Cut off misaligned data at the end
  veri_11 = [ zeros(shifting_vector(2,1)*-1,1); veri_11(1:end-(shifting_vector(2,1)*-1))];
elseif (shifting_vector(2,1) > 0)
    % OB should shift to the left, thus appending zeros at the end
    % Cut off misaligned data at the start
    veri_11 = [veri_11(1:shifting_vector(1,1)); veri_11(shifting_vector(1,1)+shifting_vector(2,1):end)];
    veri_11(numel(veri_22)) = 0;
elseif  (shifting_vector(2,1) == 0)
end
  
if size(shifting_vector,2) > 1
n = length(shifting_vector);
for i = 2:n 
    element = shifting_vector(1,i);
    delay = shifting_vector(2,i);  

    if ~(delay == 0)
        veri_11 = [veri_11(1:element-1); veri_11(element-1)*ones(delay*-1,1); veri_11(element:end)];
        veri_22(numel(veri_11)) = 0;
    end
    
end
end
close all
figure; set(gcf, 'Position', [800, 400, 700, 350]); hold on; plot(ds1); plot(ds2); hold off;  title('show shifted signals 1 and 2 for state PHI');
figure; set(gcf, 'Position', [800, 50, 700, 350]); hold on; plot(veri_11); plot(veri_22); hold off;  title('show veri signals 1 and 2 for state THETA');
end
