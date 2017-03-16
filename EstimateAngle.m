function angle = EstimateAngle(P,JPEG)
% -------------------------------------------------------------------------
%   Written March 16, 2017 by Chenglong Chen
%   Copyright 2017 by Chenglong Chen
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The end-user understands that the program was
% developed for research purposes and is advised not to rely exclusively on
% the program for any reason.
% -------------------------------------------------------------------------
% Contact: c.chenglong@gmail.com
% -------------------------------------------------------------------------
% Input:       P .... input image spectrum, e.g., Txx
%           JPEG .... flag for whether the image undergone JPEG compression
% Output:  angle .... estimated angle vector of 4-D
%                     - angle(1): using the proposed method
%                     - angle(4): using Padin et al.'s method
% -------------------------------------------------------------------------
% Estimated rotation angle based on peak found in the image spectrum using
% various methods.
% [1]
% [2] 
% -------------------------------------------------------------------------


[~,N]=size(P);

% avoiding compoments near DC with square window of size 9 x 9 (8 x 8 indeed)
% shift it to [-0.5,0.5]^2
P = fftshift(P);
% remove components near DC
w = 4;
x0 = N/2;
y0 = N/2;
P(y0-w+1:y0+w,x0-w+1:x0+w)=0;
% shift it back to [0,1]^2
P = fftshift(P);
% --------------------------------------

angle = zeros(1,4);

% nearest neighbor around the arcs in the sampling grid
Neighbor = 4;

% estimate the angle with frequency normalized to [0,1)^2
angle(1:2) = ArcEstimate_1(P,JPEG,Neighbor);
% estimate the angle with frequency normalized to [-0.5,0.5)^2
angle(3) = ArcEstimate_05(P,JPEG,Neighbor);
% estimate the angle using Padin et al.'s global searching
angle(4) = GlobalEstimate(P,JPEG);

% transform to positive angle
angle = angle+90;
end



% ----------------------------------------------------------------------- %
function EstimatedAngle = ArcEstimate_1(P,JPEG,Neighbor)
    [~,N]=size(P);
    
%     % using the symmetric of DFT
%     for k = 1
%         P = rot90(P,-k);
%     end
    % peak searching around arcs  -----------------------------------------
    theta_step = 0.2;
    AngleRange = -89:theta_step:-1;
    c = cosd(AngleRange)';
    s = sind(AngleRange)';
    
    NormalizedFlag = 1;

    % -- arc 1: l^(1,0)
    x0 = c;
    y0 = -s; 
    [X1,Y1] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % -- arc 2: l^(-1,0)
    x0 = 1-c;
    y0 = 1+s;
    [X2,Y2] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % -- arc 3: l^(0,1)
    x0 = 1+s;
    y0 = c;
    [X3,Y3] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % -- arc 4: l^(0,-1)
    x0 = -s;
    y0 = 1-c;
    [X4,Y4] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % nearest neighbor of arcs in the sampling grid  ----------
    Fx = [X1;X2;X3;X4];
    Fy = [Y1;Y2;Y3;Y4];
    
    % remove duplicate rows  ----------
    Fxy = unique([Fx,Fy], 'rows');
    Fx = Fxy(:,1);
    Fy = Fxy(:,2);
    
    % get the magnitude aroud the arcs  ----------
    linearInd = sub2ind([N,N],Fy,Fx);
    Q = P(linearInd);
    
    % peak searching  ----------
    [fx,fy] = FindPeak(Q,N,Fx,Fy,JPEG);
    
    % estimate the angle  ----------
    % find the origin of the circle
    rad = zeros(1,4);
    % (0,0)
    rad(1) = sqrt(fx^2+fy^2);
    % (1,1)
    rad(2) = sqrt((1-fx)^2+(1-fy)^2);
    % (1,0)
    rad(3) = sqrt((1-fx)^2+fy^2);
    % (0,1)
    rad(4) = sqrt(fx^2+(1-fy)^2);
    err = abs(rad-1);
    [~,ind] = min(err);
    if(ind==1)
        f1 = fx;
        f2 = fy;
    elseif(ind==2)
        f1 = 1-fx;
        f2 = 1-fy;
    elseif(ind==3)
        f1 = fy;
        f2 = 1-fx;
    else
        f1 = 1-fy;
        f2 = fx;
    end
    EstimatedAngle(1) = -atand(f2/f1);
    if(abs(EstimatedAngle(1))<=45)
        EstimatedAngle(2) = -asind(f2);
    else
        EstimatedAngle(2) = -acosd(f1);
    end
   
end



% ----------------------------------------------------------------------- %
function EstimatedAngle = ArcEstimate_05(P,JPEG,Neighbor)
    [~,N]=size(P);
    
    
%     % using the symmetric of DFT
%     for k = 1
%         P = P + rot90(P,-k);
%     end
    % shift it to [-0.5,0.5)^2
    P = fftshift(P);

    theta_step = 0.2;
    AngleRange = -89:theta_step:-1;
    c = cosd(AngleRange)';
    s = sind(AngleRange)';
    
    NormalizedFlag = 0.5;
    
    
    % -- arc 1: l^(1,0)
    x0 = c-round(c);
    y0 = -s-round(-s);
    [X1,Y1] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % -- arc 2: l^(-1,0)
    x0 = -c-round(-c);
    y0 = s-round(s);
    [X2,Y2] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % -- arc 3: l^(0,1)
    x0 = s-round(s);
    y0 = c-round(c);
    [X3,Y3] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    % -- arc 4: l^(0,-1)
    x0 = -s-round(-s);
    y0 = -c-round(-c);
    [X4,Y4] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag);
    
    
    % nearest neighbor of arcs in the sampling grid  ----------
    Fx = [X1;X2;X3;X4];
    Fy = [Y1;Y2;Y3;Y4];
    
    % remove duplicate rows  ----------
    Fxy = unique([Fx,Fy], 'rows');
    Fx = Fxy(:,1);
    Fy = Fxy(:,2);
    
    % get the magnitude aroud the arcs  ----------
    linearInd = sub2ind([N,N],Fy,Fx);
    Q = P(linearInd);
    
    % peak searching  ----------
    [fx,fy] = FindPeak(Q,N,Fx,Fy,JPEG);
    
    alpha1 = fx-0.5;
    alpha2 = fy-0.5;
    
    % [alpha1,alpha2]
    beta=atan(alpha2/alpha1);
    beta=mod(beta,pi/2);
    
    if(beta>=0&&beta<=pi/12)
        theta=-2*beta;
    elseif(beta>pi/12&&beta<=5*pi/12)
        k=cos(beta)^2*(sqrt(2*tan(beta))-tan(beta)+tan(beta)^2);
        theta=-acos(k);
    else
        theta=pi/2-2*beta;
    end
    
    EstimatedAngle = theta*180/pi;
end


% ----------------------------------------------------------------------- %
function [X,Y] = Compute_Arc_Neighbor(x0,y0,N,Neighbor,NormalizedFlag)

if(NormalizedFlag==0.5)
    x0 = x0+0.5;
    y0 = y0+0.5;
end
    if(Neighbor==2)
    % 2 nearest neighbor of x0 in the sampling grid
    x1 = ceil(x0*N);
    x2 = floor(x0*N);
    
    % 2 nearest neighbor of y0 in the sampling grid
    y1 = ceil(y0*N);
    y2 = floor(y0*N);
    
    % avoiding out of boundry
    x1(x1<1) = 1;
    x2(x2<1) = 1;
    y1(y1<1) = 1;
    y2(y2<1) = 1;
    
    % 2x2 nearest neighbor of arc l^(1,0) in the sampling grid
    X = [x1;x1;x2;x2];
    Y = [y1;y2;y1;y2];
    elseif(Neighbor==3)
            % 3 nearest neighbor of x0 in the sampling grid
    x1 = round(x0*N);
    x2 = x1-1;
    x3 = x1+1;
    
    % 3 nearest neighbor of y0 in the sampling grid
    y1 = round(y0*N);
    y2 = y1-1;
    y3 = y1+1;
    
    % avoiding out of boundry
    x1(x1<1) = 1;
    x2(x2<1) = 1;
    y1(y1<1) = 1;
    y2(y2<1) = 1;
    
    x3(x3>N) = N;
    y3(y3>N) = N;
    
    % 3x3 nearest neighbor of arc l^(1,0) in the sampling grid
    X = [x1;x1;x1;x2;x2;x2;x3;x3;x3];
    Y = [y1;y2;y3;y1;y2;y3;y1;y2;y3];
    elseif(Neighbor==4)
            % 4 nearest neighbor of x0 in the sampling grid
    x1 = ceil(x0*N)+1;
    x2 = x1-1;
    x3 = floor(x0*N);
    x4 = x3-1;
    
    % 4 nearest neighbor of y0 in the sampling grid
    y1 = ceil(y0*N)+1;
    y2 = y1-1;
    y3 = floor(y0*N);
    y4 = y3-1;
    
    % avoiding out of boundry
    x1(x1>N) = N;
    x2(x2<1) = 1;
    x3(x3<1) = 1;
    x4(x4<1) = 1;
    
    y1(y1>N) = N;
    y2(y2<1) = 1;
    y3(y3<1) = 1;
    y4(y4<1) = 1;
    
    % 4x4 nearest neighbor of arc l^(1,0) in the sampling grid
    X = [x1;x1;x1;x1;x2;x2;x2;x2;x3;x3;x3;x3;x4;x4;x4;x4];
    Y = [y1;y2;y3;y4;y1;y2;y3;y4;y1;y2;y3;y4;y1;y2;y3;y4];
    end
end


% ----------------------------------------------------------------------- %
function EstimatedAngle = GlobalEstimate(P,JPEG)
    
    [~,N]=size(P);

    % peak searching  ----------
    Q = fftshift(P);
    axes_f= 1:N;
    [Fx,Fy]=meshgrid(axes_f);
    [fx,fy] = FindPeak(Q,N,Fx,Fy,JPEG);
    
    % [alpha1,alpha2]
    alpha1 = fx-0.5;
    alpha2 = fy-0.5;
    beta=atan(alpha2/alpha1);
    beta=mod(beta,pi/2);
    
    if(beta>=0&&beta<=pi/12)
        theta=-2*beta;
    elseif(beta>pi/12&&beta<=5*pi/12)
        k=cos(beta)^2*(sqrt(2*tan(beta))-tan(beta)+tan(beta)^2);
        theta=-acos(k);
    else
        theta=pi/2-2*beta;
    end
    
    EstimatedAngle = theta*180/pi;
end


function [fx,fy] = FindPeak(Q,N,X,Y,JPEG)
    
    [~,ind]=sort(Q(:),'descend');
    if(JPEG==0)
        
        fx = X(ind(1));
        fy = Y(ind(1));
        
        % normalizing to [0,1)^2
        fx = (fx-1)/N;
        fy = (fy-1)/N;
    else
        
        L=length(Q);        
        for i=1:L
            
            fx = X(ind(i));
            fy = Y(ind(i));
            
            % normalizing to [0,1)^2
            fx = (fx-1)/N;
            fy = (fy-1)/N;
            
            fxy = [fx,fy];
            
            % ingore JPEG peaks located around (k/8,l/8)
            temp_x=round(fx*8)/8;
            temp_y=round(fy*8)/8;
            temp=[temp_x,temp_y];
            threshold=2*sqrt((1/N)^2+(1/N)^2);
            threshold=sqrt((1/N)^2+(1/N)^2);
%             threshold = 1/N;
            if(norm(temp-fxy,2)>threshold)
                break;
            end
        end
    end
end