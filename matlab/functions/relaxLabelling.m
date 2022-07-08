function [ drho ] = relaxLabelling( dX, dY, dSNR, options )
%RELAXLABELLING 
if nargin<4
    options.iterMax = 5e2;
    options.nomatch = false;
    options.convFac = 1;        % convergence power of the iteration
    options.mag = false;        % include magnitude in confidence function
end

%% make sure no NaN's are included
dX(isnan(dX)) = 0; dY(isnan(dY)) = 0; dSNR(isnan(dSNR)) = 0;

%% reweight, being the translation from SNR/corrcoef to probabilities

% include no-match label
if options.nomatch
    dX = cat(3, dX, 0.*dX(:,:,1));
    dY = cat(3, dY, 0.*dY(:,:,1));

    % no-match label, the minimum of all alternatives
    dSNR = cat(3, dSNR, nanmin(dSNR, [], 3));
end
% drho = dSNR./repmat(nansum(dSNR,3),1,1,size(dSNR,3));
drho = dSNR./repmat(sum(dSNR,3,'omitnan'),1,1,size(dSNR,3));

%% look at spatialconsistancy time

% estimate the sum of the heighest score within the neighborhood
steps = [0 0 0; ...                     % itself
    0 1 0; 0 -1 0; 1 0 0; -1 0 0; ...   % 4 direct neighbors
    1 1 0; 1 -1 0; -1 -1 0; -1 1 0];    % 4 other loose neighbors

sumCij = [];
for i = 1:size(steps, 1) % spatial neighborhood       
    cij = circshift(dSNR(:,:,1), steps(i,:));

    if steps(i,2)==1        % shift right
        cij(:,1) = NaN;
    elseif steps(i,2)==-1   % shift left
        cij(:,end) = NaN;
    end

    if steps(i,1)==1        % shift down
        cij(1,:) = NaN;
    elseif steps(i,1)==-1   % shift up
        cij(end,:) = NaN;
    end
%     sumCij = nansum(cat(3, sumCij, cij),3);
    sumCij = sum(cat(3, sumCij, cij),3,'omitnan');
end

%% main - updating of weights


for i = 1:options.iterMax
    [drhoNew] = vecCompat(dX,dY,drho,dSNR, sumCij, ...
        options.convFac, options.mag);
    drho = drhoNew; % update
end
