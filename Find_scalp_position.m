clc;
%%%%%%%%%%%%%% Instructions %%%%%%%%%%%%%%%%%%
% 1) Have the mri in workspace as mri
% 2) Have the head mask in workspace as head
% 3) Get sgACC Points
% 4) Have nasion and inion coordinates
% 5) Select chop y & z
% 6) zmid (point to get saggital circulature figure(4.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nasion = mri.SCS.NAS;
chopy = 150;
chopz = 50;


% % Colin
% nasion = [92 208.6 30]; % XYZ
inion = [91 8 50]; % XYZ
plt = [87 167 75]; plb = [87 160 67]; % sgACC left/right/top/bottom ps
prt = [95 167 73]; prb = [95 158 66];
% chopy = 150;
% chopz = 50;
% zmid = 65;


%% Getting Scalp Nodes and their Normals
h_vert = head.Vertices;
h_faces = head.Faces;

%Getting Normals
h_norm  = vert2normals(h_vert,h_faces)*-1;
% h_norm = head.VertNormals; % Does not work in regular mri

% convert MNI to MRI, need to have mri anatomical exported to workspace as "mri"
P_mri = cs_convert(mri, 'scs', 'mri', h_vert)*1000; % Scalp Nodes is voxel MRI space
N_mri = cs_convert(mri, 'scs', 'mri', h_norm)*1000; % Head Node Normals in voxel MRI space

% If the normals look weird, use this as a second method
% N_mri  = vert2normals(P_mri,h_faces)*-1; % Previous line works better.


%% Remove undersired Scalp Nodes (easier visulization)
valid_indices = P_mri(:,2) >= chopy & P_mri(:,3) >= chopz; 
h_vert = P_mri(valid_indices,:); % Scalp Nodes of Interest
h_norm = N_mri(valid_indices,:); % Scalp Nodes Normals of Interest

norming = sqrt(sum(h_norm.^2, 2));
h_norm = h_norm./ repmat(norming,1,3);

points = h_vert;
normals = h_norm;

% Show Nodes of Interest
figure(1); clf;
hold on;
% scatter3(points(:,1), points(:,2), points(:,3), 'filled', 'MarkerFaceColor', 'b');

h1 = patch('Vertices', P_mri, 'Faces', h_faces,'FaceColor', 'cyan', 'EdgeColor', 'black');
h1.FaceAlpha = .2; h1.EdgeAlpha = .1;

quiver3(points(:,1), points(:,2), points(:,3), ...
    normals(:,1), normals(:,2), normals(:,3), 'r');

xlabel('X'); ylabel('Y'); zlabel('Z'); view(3); axis equal;


%% Getting brain ROI lines 
%points of interest
p11 = plt; p12 = plb;
p21 = prt; p22 = prb;

line1 = get_3d_line(plt,plb);
line2 = get_3d_line(prt,prb);


%% Getting Lines Projections of scalp Nodes
% This makes the projections of each Scalp Node into the head.

distance = 66; % projection distance (mm)
density = 2; % number of points per 1mm
N = distance*2+1; % Number of points for each projection, This indicates how deep the projections go
z = linspace(0,distance,N).';
Nl = length(h_vert);
lines = zeros(N,3,Nl);

for i = 1:Nl
    p = h_vert(i,:);
    dir = normals(i,:)*(-1);
    lines(:,:,i) = repmat(p,N,1) + z.*dir; 
end


%% Getting Colsest and Best Scalp Points;
% Not only closest, but also one that best sonicates the area of interest.

errors = zeros(Nl,2);
l1 = length(line1); l2 = length(line2);

for i = 1:Nl
    
    p = lines(:,:,i);
    err = zeros(l1,1);
    for i1 = 1:l1
        d1 = repmat(line1(i1,:),N,1);
        tmp = p-d1;
        tmp = sqrt(sum(tmp.^2,2));
        err(i1) = min(tmp);
    end
%     errors(i,1) = min(abs(err)); %min distance
    errors(i,1) = sum(err); % min distance with best aligment.
    err = zeros(l2,1);
    for i2 = 1:l2
        d2 = repmat(line2(i2,:),N,1);
        tmp = p-d2;
        tmp = sqrt(sum(tmp.^2,2));
        err(i2) = min(tmp);
    end
%     errors(i,2) = min(abs(err)); %min distance
    errors(i,2) = sum(err); % min distance with best aligment.
end
        
[min1,I1] = min(errors(:,1)); % Find min sum error for ROI 1
[min2,I2] = min(errors(:,2)); % Find min sum error for ROI 2
        
p1 = lines(1,:,I1); % Find Best Scalp Point for ROI 1
p2 = lines(1,:,I2); % Find Best Scalp Point for ROI 2

tmp = repmat(p1,l1,1) - line1;
distance_1 = sqrt(sum(tmp.^2,2));
distance_1 = min(distance_1);

tmp = repmat(p1,l2,1) - line2;
distance_2 = sqrt(tmp(:,1).^2+tmp(:,2).^2+tmp(:,3).^2);
distance_2 = min(distance_2);

fprintf('\nP1_left_sgACC = %.3f mm',distance_1);
fprintf('\nP2_right_sgACC = %.3f mm',distance_2);
        
    
%% Plotting Best Scalp Points

figure(2); clf; hold on;
points = permute(lines,[1 3 2]);
points = reshape(points,[],3);
 
plot3(points(:,1), points(:,2), points(:,3),'Color', [0.2 0 0 0.1]);
view(3);

points = line1;
plot3(points(:,1), points(:,2), points(:,3),'b*'); axis equal;
points = line2;
plot3(points(:,1), points(:,2), points(:,3),'m*');
xlabel('X'); ylabel('Y'); zlabel('Z'); view(3);

points = lines(:,:,I1);
plot3(points(:,1), points(:,2), points(:,3),'r.'); axis equal;
points = lines(:,:,I2);
plot3(points(:,1), points(:,2), points(:,3),'g.');

plot3(p1(1),p1(2),p1(3),'.r','MarkerSize',20);
plot3(p2(1),p2(2),p2(3),'.g','MarkerSize',20);


%% Ploting Head and Best Scalp Points
    
figure(3); clf; axis equal; hold on;
points = P_mri;

h1 = patch('Vertices', points, 'Faces', h_faces, ...
      'FaceColor', 'cyan', 'EdgeColor', 'black');
h1.FaceAlpha = .2;
h1.EdgeAlpha = .1;
  
plot3(p1(1),p1(2),p1(3),'.r','MarkerSize',20);
plot3(p2(1),p2(2),p2(3),'.g','MarkerSize',20);

points = line1;
plot3(points(:,1), points(:,2), points(:,3),'*b');
points = line2;
plot3(points(:,1), points(:,2), points(:,3),'*m');

points = lines(:,:,I1);
plot3(points(:,1), points(:,2), points(:,3),'.r'); axis equal;

points = lines(:,:,I2);
plot3(points(:,1), points(:,2), points(:,3),'.g'); axis equal;

xlabel('X'); ylabel('Y'); zlabel('Z');

view(3);
hold off;    
    

%% Measuring Curvature.

% values for circsort
xmid = nasion(1);
ymid = (max(P_mri(:,2))+min(P_mri(:,2)))/2;
zmid = (inion(:,3)+nasion(3))/2;

xsp = 1;
zsp = 1;

% Getting  Mid Saggital Nodes
points = P_mri; % Nodes points of the MRI
valid_indices = points(:,1) >= xmid-xsp & points(:,1) <= xmid+xsp; 
points = points(valid_indices,:);
points = points(:,2:3);

% Sorting nodes (allows curvature measurement)
points3 = circsort(points,[xmid zmid],+90,1); 

% project point P in the scalp curvature and add it to the list
[points3,pp1] = coronal_projection(points3,p1(2:3));
[points3,pp2] = coronal_projection(points3,p2(2:3));
[points3,pn] = pjp2cont(points3,nasion(2:3));
[points3,pi] = pjp2cont(points3,inion(2:3));

points3 = circsort(points3,[xmid zmid],+90,0); %angle sort
dis1 = circlen(points3,pp1,pn); % point1 Distance
dis2 = circlen(points3,pp2,pn); % point2 Distance
circ1 = circlen(points3,pi,pn);

fprintf('\n\nNasion2Inion = %.1f mm',circ1);
fprintf('\nP1, P2 = (%.1f, %.1f) mm',dis1,dis2);
fprintf('\nP1, P2 = (%.1f%%, %.1f%%) ',dis1/circ1*100,dis2/circ1*100);

figure(4); clf;
subplot(131); hold on; axis equal; axis([-10 230 -10 180]);
plot(points3(:, 1), points3(:, 2), 'k');
plot(nasion(2),nasion(3),'.b','MarkerSize',20)
plot(inion(2),inion(3),'.m','MarkerSize',20)
plot(pp1(1),pp1(2),'*r','MarkerSize',20);
plot(p1(2),p1(3),'.r','MarkerSize',20);
plot(pp2(1),pp2(2),'*g','MarkerSize',20);
plot(p2(2),p2(3),'.g','MarkerSize',20);
title('saggital view'); legend('','Nasion','Inion','P1 Projection','P1','P2 Projection','P2')

% For p1
z_loc = p1(3);
points = P_mri; 
valid_indices = points(:,3)>=z_loc-zsp & points(:,3)<=z_loc+zsp;
points = points(valid_indices,:);
points = points(:, 1:2);
points3 = circsort(points,[xmid ymid],-90,1); %angle sort
% project and add points
[points3,pp1] = pjp2cont(points3,p1(1:2));
[points3,pn] = pjp2cont(points3,nasion(1:2));
points3 = circsort(points3,[xmid ymid],-90,0); %angle sort
dis = circlen(points3,pp1,pn); % point Distance 
points3 = [points3(:,1:2); points3(1,1:2)];
circ1 = circlen2(points3);

fprintf('\n\nAxial Circ = %.1f mm',circ1);
fprintf('\nP1 = %.1f mm',dis);
fprintf('\nP1 %.1f%%',dis/circ1*100);

subplot(132); hold on; axis equal; axis([-10 200 0 220]);
plot(points3(:,1), points3(:,2), 'k');
plot(pn(1), pn(2),'b*','MarkerSize',20);
plot(nasion(1), nasion(2),'.b','MarkerSize',20);
plot(pp1(1),pp1(2),'*r','MarkerSize',20)
plot(p1(1),p1(2),'.r','MarkerSize',20);
set(gca, 'YDir','reverse', 'XDir','reverse')
title('Axial P1'); legend('','Nasion Projection','Nasion','P1 Projection','P1')

% For p2
z_loc = p2(3);
points = P_mri; 
valid_indices = points(:,3)>=z_loc-zsp & points(:,3)<=z_loc+zsp;
points = points(valid_indices,:);
points = points(:, 1:2);
points3 = circsort(points,[xmid ymid],-90,1); %angle sort
% project and add points
[points3,pp2] = pjp2cont(points3,p2(1:2));
[points3,pn] = pjp2cont(points3,nasion(1:2));

points3 = circsort(points3,[xmid ymid],-90,0); %angle sort
% Point Distance
dis = circlen(points3,pp2,pn); 
points3 = [points3(:,1:2); points3(1,1:2)];
circ1 = circlen2(points3);

fprintf('\n\nAxial Circ = %.1f mm',circ1);
fprintf('\nP2 = %.1f mm',dis);
fprintf('\nP2 %.1f%%\n',dis/circ1*100);

subplot(133); hold on; axis equal; axis([-10 200 0 220]);
plot(points3(:,1), points3(:,2), 'k');
plot(pn(1), pn(2),'b*','MarkerSize',20);
plot(nasion(1), nasion(2),'b.','MarkerSize',20);
plot(pp2(1),pp2(2),'*g','MarkerSize',20)
plot(p2(1),p2(2),'.g','MarkerSize',20);
set(gca, 'YDir','reverse', 'XDir','reverse')
title('Axial P2'); legend('','Nasion Projection','Nasion','P2 Projection','P2')




%% functions

function [ordXY] = circsort(points,center,angle_delay,smoothe)
% Sort unorgonized 2d points by their angle to a center point
    
    vectors = points - repmat(center,length(points),1);
    angles = atan2(vectors(:,2),vectors(:,1))*180/pi+angle_delay;
    angles = mod(angles, 360);
    points(:,3) = angles;
    sorted_p = sortrows(points, 3);
    
    ordXY = sorted_p(:,1:2);
    
    if smoothe == 1 % Smoothing with Moving Average
        windowSize = 5; % Size of the moving window
        smoothX = movmean(ordXY(:,1),windowSize);
        smoothY = movmean(ordXY(:,2),windowSize);
        ordXY = [smoothX, smoothY];
    end
end


function dis = circlen(coor,p,ref)
% Find curvature between point of interest and nasion.
    
    idp = find(coor == p); idp = idp(1);
    idref = find(coor == ref); idref = idref(1);
    
    a = min([idp, idref]); b = max([idp, idref]);
    c = length(coor)+a;
    
    cor2 = coor(a+1:b,:); cor1 = coor(a:b-1,:);
    tmp = cor2-cor1;
    tmp = sqrt(sum(tmp.^2,2));
    dis = sum(tmp);
 
    
    if coor(end,:) ~= coor(1,:)
        coor = [coor; coor];
    else
        coor = [coor; coor(2:end,:)];
    end
    
    cor2 = coor(b+1:c,:); cor1 = coor(b:c-1,:);
    tmp = cor2-cor1;
    tmp = sqrt(sum(tmp.^2,2));
    dis2 = sum(tmp);

    if dis2 < dis
        dis = dis2;
    end
    
end

function dis = circlen2(coor)
% Find Axial Curvature
    if coor(end,:) ~= coor(1,:)
        coor(end+1) = coor(1,:);
    end
    
    cor2 = coor(1+1:end,:); cor1 = coor(1:end-1,:);
    tmp = cor2-cor1;
    tmp = sqrt(sum(tmp.^2,2));
    dis = sum(tmp);
end

function [contour_plus,pp] = pjp2cont(coordinates,p)
% Project Point p into 2d Contour

    % New version, Vectorized to N points
    N = 7;
    tmp = coordinates - repmat(p, length(coordinates),1);
    tmp = sqrt(sum(tmp.^2,2));
    
    mins = zeros(N,1);
    for i = 1:N
        [~,mins(i)] = min(tmp); tmp(mins(i)) = inf;
    end
    
    mins_coors = coordinates(mins,:);
    A = mins_coors(1,:);
    mins_coors_test = mins_coors(2:end,:);
    
    AP = p - A; AP = repmat(AP,N-1,1);
    points_dis = mins_coors_test - repmat(A,N-1,1);
    
    t_vars = dot(AP,points_dis,2)./dot(points_dis,points_dis,2);
%     t_vars(t_vars>1) = 1; t_vars(t_vars<0) = 0;
    
    proj_vars = A + t_vars .* points_dis;
    
    dist_projections = norm(repmat(p,N-1,1) - proj_vars);
    
    [~,id] = min(dist_projections);
    
    pp = proj_vars(id,:);
    
%     tmp = coordinates - repmat(p, length(coordinates),1);
%     tmp = sqrt(sum(tmp.^2,2));
%     
%     [~,min1] = min(tmp); tmp(min1) = inf;
%     [~,min2] = min(tmp); tmp(min2) = inf;
%     [~,min3] = min(tmp);
%     
%     A = coordinates(min1,:); B = coordinates(min2,:);
%     C = coordinates(min3,:);
%     
%     AP = p - A; AB = B - A; AC = C - A;
%     
%     tB = dot(AP, AB)/dot(AB, AB); 
%     tB(tB>1) = 1; tB(tB<0) = 0;
%     projB = A + tB * AB;
%     
%     tC = dot(AP, AC)/dot(AC, AC);
%     tC(tC>1) = 1; tC(tC<0) = 0;
%     projC = A + tC * AC;
% 
%     dist_AB = norm(p - projB);
%     dist_AC = norm(p - projC);
% 
%     if dist_AB <= dist_AC
%         pp = projB;
%     elseif dist_AC <= dist_AB
%         pp = projC;
%     else
%         % pass
%     end
    
    contour_plus = [pp; coordinates];  
    contour_plus = unique(contour_plus,'rows');   
end

function [contour_plus,pp] = coronal_projection(coordinates,p)
% Project Point p coronal axis into 2d Contour

    tmp = abs(coordinates(:,2) - p(2));
    
    c = zeros(4,1); tr = []; tl = [];
    fail_safe = 0;
    keep_search = 1;
    while keep_search
       
       fail_safe = fail_safe + 1;
       if fail_safe > length(coordinates)*2
           break
       end
       
       [~, id] = min(tmp);
%        t = [t; coordinates(id,:)];
       tmp(id) = inf;
       
       for i = 1:length(id)
           
           ii = id(i);
           t = coordinates(ii,:);
           if t(1) >= median(coordinates(:,1))
               if t(2) >= p(2) && c(1) == 0
                   tr = [tr; t]; c(1) = 1;
               elseif t(2) < p(2) && c(2) == 0
                   tr = [tr; t]; c(2) = 1; 
               end
           elseif t(1) < median(coordinates(:,1))
               if t(2) >= p(2) && c(3) == 0
                   tl = [tl; t]; c(3) = 1;
               elseif t(2) < p(2) && c(4) == 0
                   tl = [tl; t]; c(4) = 1; 
               end
           end
                   
        if sum(c) == 4
            keep_search = 0; end
       end
    end
    
    % frontal projection
    if size(tr,1) == 2
        m = (tr(2,2) - tr(1,2))/(tr(2,1) - tr(1,1));
        b = tr(2,2) - m*tr(2,1);
        x = (p(2) - b) / m;
        p_anterior = [x, p(2)];
    end
    
    % frontal projection
    if size(tl,1) == 2
        m = (tl(2,2) - tl(1,2))/(tl(2,1) - tl(1,1));
        b = tl(2,2) - m*tl(2,1);
        x = (p(2) - b) / m;
        p_posterior = [x, p(2)];
    end

    if exist('p_anterior') && exist('p_posterior')
        if norm(p-p_anterior) <= norm(p-p_posterior)
            pp = p_anterior;
        else
            pp = p_posterior;
        end
    elseif exist('p_anterior')
        pp = p_anterior;
    elseif exist('p_posterior')
        pp = p_posterior;
    end
        
    
    contour_plus = [pp; coordinates];  
    contour_plus = unique(contour_plus,'rows');   
end

function normals = vert2normals(vertices,faces)
    % Calculate the normals for each face and accumulate them at the corresponding vertices

    vertex_normals = zeros(length(vertices), 3);

    for i = 1:length(faces)
        % Get the vertex indices of the i-th face
        face = faces(i, :);

        % Extract the vertices of the current face
        v1 = vertices(face(1), :);
        v2 = vertices(face(2), :);
        v3 = vertices(face(3), :);

        % Calculate the normal of the current face using cross product
        normal = cross(v2 - v1, v3 - v1);

        % Accumulate the normal at each vertex of the face
        vertex_normals(face, :) = vertex_normals(face, :) + normal;
    end

    % Normalize vertex normals
    normals = vertex_normals ./ sqrt(sum(vertex_normals.^2, 2));

    % 'vertex_normals' now contains the normal vectors at each vertex

end

function line = get_3d_line(p1,p2)

    N = ceil(norm(p1-p2));
    x = linspace(p1(1), p2(1), N);
    y = linspace(p1(2), p2(2), N);
    z = linspace(p1(3), p2(3), N); 
    line = [x; y; z]';
    
    if p1 == p2
        line = p1;
    end
end