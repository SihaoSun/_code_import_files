function Data_filtered = bebop_prefilter(DataStruct)

DStmp = DataStruct;
fields = fieldnames(DStmp);

delta_T = DStmp.TIME(2)-DStmp.TIME(1);
Ndata = length(DStmp.TIME);


% Remove outliers from dataset by investigating the behavior of the dataset derivatives
OutlierCount = 0;
CorrectedCount = 0;
CorrectedFieldCount = cell(length(fields), 2);
stdBounds = 20; % maximum number of allowed standard deviations in the signal
stdBoundsD = 1; % maximum number of allowed standard deviations in the derivative of the signal
searchWnd = 6;

    for j = 1:length(fields)-1
        field   = fields{j};
        if ((strcmp(field,'att_G2O') || strcmp(field,'posCO_G') || strcmp(field,'posCG_E') ||...
              strcmp(field,'vel_E') || strcmp(field,'vel_B') || strcmp(field,'velCG_B') ||...
              strcmp(field,'accCG_E') || strcmp(field,'accCG_B') || strcmp(field,'OMEGA_B') ||...
              strcmp(field,'angacc_B') || strcmp(field,'velCG_E') || strcmp(field,'acc_E') ||...
              strcmp(field,'acc_B') || strcmp(field,'AoA') || strcmp(field,'AoA_b') ||...
              strcmp(field,'beta_b') || strcmp(field,'ay_b') || strcmp(field,'az_b') ||...
              strcmp(field,'TAS') || strcmp(field,'IAS') || strcmp(field,'Mach') || strcmp(field,'tat') ||...
              strcmp(field,'palt') || strcmp(field,'bcalt') || strcmp(field,'thd')))
            continue;
        end



        CorrectedFieldCount{j, 1} = field;
        CorrectedFieldCount{j, 2} = 0;
        FixedCount{j, 1}     = field;
        FixedCount{j, 2}     = 0;
        if (isempty(DStmp.(field)))
            continue;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % plot data
%         figure;
%         hold on;
%         plot(DStmp.TIME, DStmp.(field));
%         ylabel(field);
%         xlabel('time [s]');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % get statistics

        while 1
            DStmp.(field) = interpNan(DStmp.(field));
            meanf   = mean(DStmp.(field));
            stdf    = std(DStmp.(field));
            % find outliers by searching for any points in the dataset itself that are more than 'stdBounds' standard deviations from the mean
            outliersFn = find(abs(DStmp.(field)-meanf) > stdBounds*stdf);
            if ~isempty(outliersFn)
                for k = length(outliersFn)
                    i = 1;
                    while ~isempty(intersect(outliersFn(k)-i,outliersFn))
                        i = i + 1;
                    end
                    DStmp.(field)(outliersFn(k)) = DStmp.(field)(outliersFn(k)-i);
                end
            else
                break;
            end
        end
        df      = diff(DStmp.(field));
        stddf   = std(df);
        meandf  = mean(df);        
        % find outliers by searching for any points in the derivative of the dataset that are more than 'stdBoundsD' standard deviations from the mean
        outliers = find(abs(df-meandf) > stdBoundsD*stddf);
        outliers = unique([outliers; outliersFn]);
        if (~isempty(outliers))
            % go through all outliers, and correct them
            OutlierCount = OutlierCount + length(outliers);
            for k = 1:length(outliers)
                ol = outliers(k);
                if (ol > searchWnd)
                    % find candidate replacement within 'searchWnd' samples from the outlier
                    candfound = 0;
                    bestcand = ol-1; % Estimate of the best candidate
                    wndupper = 3;
                    if (ol+wndupper > Ndata)
                        wndupper = Ndata-ol;
                    end
                    maxdfold = max(abs(diff( DStmp.(field)(ol-wndupper:ol+wndupper) )));
                    min_maxdftmp = inf;
                    for olc = ol-1:-1:ol-searchWnd
                        if (abs(df(olc)-meandf) < stdBoundsD*stddf)

                            % test the candidate: the resulting derivative should be smaller than first
                            candval = DStmp.(field)(olc+1);
                            tmp = DStmp.(field);
                            tmp(ol+1) = candval;
                            dftmp = diff(tmp(ol-wndupper:ol+wndupper));
                            maxdftmp = max(abs(dftmp));
                            if (maxdftmp <= min_maxdftmp)
                                min_maxdftmp = maxdftmp;
                                bestcand = olc;
                                candfound = 1;
                            end

                            %break;
                        end
                    end
                    if (candfound)

                        candval = DStmp.(field)(bestcand+1);
                         % the new candidate reduces the maximum derivative!
                        DStmp.(field)(ol+1) = candval; % + 1 because df has one sample less!
                        CorrectedFieldCount{j, 2} = CorrectedFieldCount{j, 2} + 1;
                        CorrectedCount = CorrectedCount + 1;

                        % test the candidate: the resulting derivative should be smaller than first
%                             candval = DStmp.(field)(olc+1);
%                             tmp = DStmp.(field);
%                             tmp(ol+1) = candval;
%                             wndupper = 3;
%                             if (ol+wndupper > Ndata)
%                                 wndupper = Ndata-ol;
%                             end
%                             dftmp = diff(tmp(ol-wndupper:ol+wndupper));
%                             maxdftmp = max(abs(dftmp));
%                             maxdfold = max(abs(diff( DStmp.(field)(ol-wndupper:ol+wndupper) )));
%                             
% %                            if (abs(candval-DStmp.(field)(ol)) < abs(df(ol)-meandf))
%                             if (maxdftmp <= maxdfold)
%                                 % the new candidate reduces the maximum derivative!
%                                 DStmp.(field)(ol+1) = candval; % + 1 because df has one sample less!
%                                 CorrectedFieldCount{j, 2} = CorrectedFieldCount{j, 2} + 1;
%                                 CorrectedCount = CorrectedCount + 1;
%                             end
                    end
                else
                    DStmp.(field)(ol) = DStmp.(field)(ol+1);
                    CorrectedFieldCount{j, 2} = CorrectedFieldCount{j, 2} + 1;
                    CorrectedCount = CorrectedCount + 1;
                end % if ol > 1
            end % go through all outliers

        end

%         plot(DStmp.TIME, DStmp.(field),'r--');

    end % go through all fields

    Data_filtered = DStmp;
end




