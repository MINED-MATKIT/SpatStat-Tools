function [rotstats, xx, yy] = RotInv(Cartesianstats,varargin)
%RotInv(Cartesianstats) will obtain rotationally invariant representation
%from two point statistics. Cartesianstats is a centered statistics 2D matrix.
%
%Optional:
%   RotInv(Cartesianstats,anglebindwidth) will discretize theta space
%   with the specified bin width.
%
%   RotInv(Cartesianstats,'figure') will draw a figure with the output.
%
% A by-product of ongoing computational materials science research at 
% MINED@Gatech.(http://mined.gatech.edu/)
%
% Copyright (c) 2017, Ahmet Cecen and MINED@Gatech  -  All rights reserved.

    % Initialize Optional Parameters
    anglebinwidth = 3;
    figureswitch = false;
    
    % Accept Optional Parameters
    if nargin > 1
        for ii=1:length(varargin)
            if isscalar(varargin{ii})
                anglebinwidth = varargin{ii};
            elseif strcmp(varargin{ii},'figure')
                figureswitch = true;
            else
                error('Optional argument must be a scalar or ''figure''');
            end
        end
    end

    % Find Center
    cx = floor(size(Cartesianstats,1)/2)+1;
    cy = floor(size(Cartesianstats,2)/2)+1;
    
    % Obtain Polar Coordinates
    [X, Y]=find(ones(size(Cartesianstats)));
    [theta,rho] = cart2pol(X - cx , Y - cy);

    % Find Padding Indices
    lowerpadthetaind = find(theta>=0.5*pi); 
    upperpadthetaind = find(theta<=-0.5*pi); 
    
    % Pad Polar Coordinates for Periodicity along Theta
    theta = [theta(lowerpadthetaind)-2*pi; theta; theta(upperpadthetaind)+2*pi];
    rho = [rho(lowerpadthetaind);rho;rho(upperpadthetaind)];
    gg2 = [Cartesianstats(lowerpadthetaind);Cartesianstats(:);Cartesianstats(upperpadthetaind)];

    % Accomodare Volume Fraction Information
    th3 = [theta; ([-270:-1,1:270]*pi/180)'];
    rh3 = [rho; zeros(540,1)];
    gg3 = [gg2; Cartesianstats(cx,cy)*ones(540,1)];

    % Fit Interpolant
    FFInt2 = scatteredInterpolant(th3,rh3,gg3);
    
    % Define R,Theta Grid
    [xt, yt] = meshgrid(-180:anglebinwidth:179,0:min([cx,cy]));
    xt = xt*pi/180;

    % Obtain Values on Uniform Grid
    vq2 = FFInt2(xt(:),yt(:));
    VQ2 = reshape(vq2,[min([cx,cy]+1),size(xt,2)]);

    % Impose Rotational Invariance / Zero Shift Phases
    rotstats = ifft(abs(fft(VQ2,[],2)),[],2);
    [xx, yy] = pol2cart(xt,yt);
    
    % Optional Figure
    if figureswitch
        figure('pos',[100 100 700 600]);
        surf([xx xx(:,1)],[yy yy(:,1)],[rotstats rotstats(:,1)]);
        axis image
        shading interp
        axis off
        view(2)
        colorbar
        set(gca,'FontSize',20)

        hold on

        plot3([-max(xx(:))-50 -max(xx(:))-50],[0 max(xx(:))],[1 1],'k','LineWidth',2)
        [axisx, axisy] = pol2cart((0:anglebinwidth:360)*pi/180,repelem(max(xx(:)),size(xx,2)+1));
        plot3(axisx,axisy,ones(length(axisx),1),'k','LineWidth',2);

        text(-max(xx(:))-50-5,0,1,'0','FontSize',16,'HorizontalAlignment','right')
        text(-max(xx(:))-50-5,max(xx(:)),1,int2str(max(xx(:))),'FontSize',16,'HorizontalAlignment','right')

        text(-max(xx(:))-1,0,1,'180','FontSize',16,'HorizontalAlignment','right')
        text(max(xx(:))+1,0,1,'0','FontSize',16,'HorizontalAlignment','left')
        text(0,-max(xx(:))-10,1,'270','FontSize',16,'HorizontalAlignment','center')
        text(0,max(xx(:))+13,1,'90','FontSize',16,'HorizontalAlignment','center')

        for ii = 1:floor(max(xx(:))/100)
            [axisx, axisy] = pol2cart((0:anglebinwidth:360)*pi/180,repelem(ii*100,size(xx,2)+1));
            plot3(axisx,axisy,ones(length(axisx),1),'--k','LineWidth',1);
            text(-max(xx(:))-50-5,ii*100,1,num2str(ii*100),'FontSize',16,'HorizontalAlignment','right')
            plot3([-max(xx(:))-50 -max(xx(:))-20],[ii*100-2 ii*100-2],[1 1],'--k','LineWidth',1)
        end
    end

end

