function [ dPnew ] = vecCompat( dX, dY, dP, dRho, varargin)
%VECCOMPAT estimates the new confidence of a multimodel vector estimation
% based on the knowledge of other vector fields
%lk1`.  
% dX -      vector component in X-direction
% dY -      vector component in Y-direction
% dP -      weighting per candidate
% dRho -    correlation score per candidate
% following entries are in quadruples (dX,dY,dPdRho), and
% ending with
% Cij -     matrix with weighting coeff
% w -       convergence
% mag -     boolean to include(1) or exclude(0) the influence of difference magnitude

if mod(nargin-3,4)~=0
    error('not the right number of inputs: dX, dY, dP, dRho...')
else                                    % reorder pairs into structure
    tempDim = ((nargin-3)/4)-1;
    for h = 1:tempDim
        dX2{h} = varargin{h.*4-3};
        dY2{h} = varargin{h.*4-2};
        dP2{h} = varargin{h.*4-1};
        dRho2{h} = varargin{h.*4};
    end
    sumCij = varargin{end-2};
    w = varargin{end-1};                % convergence function
    mag = varargin{end};
end

n = size(dX,3);
q = nan.*dX;

neigh = 8;
steps = [0 0 0; ...                     % itself
    0 1 0; 0 -1 0; 1 0 0; -1 0 0; ...   % 4 direct neighbors
    1 1 0; 1 -1 0; -1 -1 0; -1 1 0];    % 4 other loose neighbors

for j = 1:n
    a = cat(3, dX(:,:,j), dY(:,:,j));       % vector of interest
    totalSupport = [];
    
    %% look at vectors in the temporal domain
    for h = 1:tempDim+1                             % matchable vectors in time
        if h == 1
            opt = size(dX,3);
            rho = dRho;
        else
            opt = size(dX2{h},3);           % options to be the best
            rho = dRho2{h};
        end
        if tempDim==0
            searchSpace = 2:neigh;
        else
            searchSpace = 1:neigh;
        end
            
        for i = searchSpace;                    % matchable in spatial domain
        
            for k = 1:opt                   % alternatives
                
                if ~((h == 1) && (i == 1))
                    %% get information about neighbor
                    if h == 1               % for same timeframe
                        b = cat(3, circshift(dX(:,:,k), steps(i,:)), ...
                        circshift(dY(:,:,k), steps(i,:)));
                        phk = circshift(dP(:,:,k), steps(i,:));
                    else                    % for different timeframe
                        b = cat(3, circshift(dX2{h}(:,:,k), steps(i,:)), ...
                            circshift(dY2{h}(:,:,k), steps(i,:)));
                        phk = circshift(dP2{h}(:,:,k), steps(i,:));
                    end

                    %% boundary conditions
                  
                    if steps(h,2)==1        % shift right
                        b(:,1) = NaN; phk(:,1) = NaN;
                    elseif steps(h,2)==-1   % shift left
                        b(:,end) = NaN; phk(:,end) = NaN;
                    end

                    if steps(h,1)==1        % shift down
                        b(1,:) = NaN; phk(1,:) = NaN;
                    elseif steps(h,1)==-1   % shift up
                        b(end,:) = NaN; phk(end,:) = NaN;
                    end

                    %% estimate compatability
                    % angular function between two vector, theta = [-1 ... 1]
%                     theta = dot(a,b,3) ./ (sqrt(nansum(a.^2,3)).* sqrt(nansum(b.^2,3))); 
                    theta = dot(a,b,3) ./ (sqrt(sum(a.^2,3,'omitnan')).* sqrt(sum(b.^2,3,'omitnan'))); 

                    if mag
                        % include magnitude
%                         La = sqrt(nansum(a.^2,3)); % calculate length
%                         Lb = sqrt((nansum(b.^2,3)));
                        La = sqrt(sum(a.^2,3,'omitnan')); % calculate length
                        Lb = sqrt((sum(b.^2,3,'omitnan')));
                        L = (1 - abs(La-Lb)./max(cat(3,La,Lb),[],3)); % magnitude factor weighting

                        % estimate support
                        localSupport = theta.*L.*phk;
                    else
                        localSupport = theta.*phk;
                    end
                    localSupport(isnan(localSupport)) = 0;
                    if k == 1
                        if h == 1
                            counter = ~isnan(phk);
                        end
                        support = localSupport;
                    else
                        counter = counter + ~isnan(phk);
                        support = support+localSupport;
                    end
                end
            end
            
            %% confidence coefficient
            Cij = circshift(rho(:,:,1), steps(i,:));
            if steps(h,2)==1        % shift right
                Cij(:,1) = NaN;
            elseif steps(h,2)==-1   % shift left
                Cij(:,end) = NaN;
            end

            if steps(h,1)==1        % shift down
                Cij(1,:) = NaN;
            elseif steps(h,1)==-1   % shift up
                Cij(end,:) = NaN;
            end
            
            neighSupport = (Cij./sumCij).*support; % weight coefficients
%             totalSupport = nansum(cat(3,totalSupport, neighSupport),3);
            totalSupport = sum(cat(3,totalSupport, neighSupport),3,'omitnan');
        end
    end

    q(:,:,j) = ((1./(counter-1)) .* totalSupport); % total support
end

% dPnew = ( (dP(:,:,1:n).*(1 + q))./repmat(nansum(dP(:,:,1:n).*(1 + q), 3), 1, 1, n) ).^w;   
dPnew = ( (dP(:,:,1:n).*(1 + q))./repmat(sum(dP(:,:,1:n).*(1 + q), 3,'omitnan'), 1, 1, n) ).^w;           
end

