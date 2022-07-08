function [THETA] = relaxLabel_multiLevel_output(Theta)

primeDir = Theta(:,:,1);
seconDir = Theta(:,:,3);

dX = cat(3, -sind(2.*primeDir), -sind(2.*seconDir));
dY = cat(3, cosd(2.*primeDir), cosd(2.*seconDir));
dSNR = cat(3, Theta(:,:,2), Theta(:,:,4));

fprintf('Applying Relaxation labelling to get a smooth transition in estimate.\n');
drho = relaxLabelling(dX,dY,dSNR );

[~,idx] = max(drho,[],3);

idx(idx==2) = 3;
[M,N,B] = size(Theta);
[Nm,Mm] = meshgrid(1:N,1:M);
[f] = sub2ind([M,N,B],Mm(:),Nm(:),idx(:));
% primeDir = reshape(Theta(f),M,N)-90; % theta-90; % dominant crevasse direction (w.r.t. x-ax)
primeDir = reshape(Theta(f),M,N); % theta-90; % dominant crevasse direction (w.r.t. x-ax)
[f] = sub2ind([M,N,B],Mm(:),Nm(:),idx(:)+1);
Score1 = reshape(Theta(f),M,N); % primary crevasse signal strength [-90 90]
idx(idx==3) = 2; idx(idx==1) = 3; idx(idx==2) = 3;
[f] = sub2ind([M,N,B],Mm(:),Nm(:),idx(:));
% seconDir = reshape(Theta(f),M,N)-90; % secondary crevasse direction
seconDir = reshape(Theta(f),M,N); % secondary crevasse direction
[f] = sub2ind([M,N,B],Mm(:),Nm(:),idx(:)+1);
Score2 = reshape(Theta(f),M,N); % secondary crevasse signal strength

clear M N B f idx drho val dSNR dX dY %Theta

THETA = cat(3, primeDir, Score1, seconDir, Score2, Theta(:,:,5:end) );
% fprintf(['Estimated directions for ' imName(1:end-3) '\n']);
end
