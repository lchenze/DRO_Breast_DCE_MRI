function les_mask = LesionGenerator(obj)
% create paramaterized lesions
% 1 is round, 2 is lobulated, 3 is irregular
% and 4 is spiculated

R = obj.lesionBaseSize;

switch obj.lesionShape
    case 'round'
        type = 1;
    case 'lobulated'
        type = 2;
    case 'irregular'
        type = 3;
    case 'spiculated'
        type = 4;
end

%Set the random number generator based upon the parameters set:
if obj.randSeed == 0
    rng(round((type + obj.lesionSize + obj.lesionFeatureSize*10 + obj.lesionFeatures)*10));
else
    rng(obj.randSeed);
end

% generate base lesion points
lpts = GenerateLesion(R,type, obj);

data.x = lpts.x;
data.y = lpts.y;
data.z = lpts.z;

% get volume data
data.res = 100;
data.lesvol = ConvertToVolume(data);
les_mask = data.lesvol;

% plot resulting lesion
% LesionPlot(data)
% figure
% mon_vol = reshape(les_mask,[size(les_mask,1),...
%                 size(les_mask,2),1,size(les_mask,3)]);
% montage(mon_vol)
end

%% subfunctions %%

function lesionpoints = GenerateLesion(R,type, obj)
    % generate base lesion
    [x1, y1, z1] = sphere(150);
    x1 = R*x1(:);
    y1 = R*y1(:);
    z1 = R*z1(:);
    if type == 1
        lesionpoints = struct('x',x1,'y',y1,'z',z1);
    elseif type == 2    %lobulated
        
        % add 3 medium lobules
        numB = obj.lesionFeatures;
        theta = 2*pi*rand(1,numB);
        phi = asin(2*rand(1,numB)-1);
        [xadd,yadd,zadd] = sph2cart(theta,phi,.95*R*ones(1,numB));
        for k = 1:numB
           [x,y,z] = sphere(150);
           x = (R)*x(:) + xadd(k);
           y = (R)*y(:) + yadd(k);
           z = (R)*z(:) + zadd(k);
           x1 = [x1; x];
           y1 = [y1; y];
           z1 = [z1; z];
        end
        lesionpoints = struct('x',x1,'y',y1,'z',z1);
        
%         %add 3 lobules
%         xshift = .8*R*(rand(obj.lesionFeatures,1))-.4*R;
%         yshift = .8*R*(rand(obj.lesionFeatures,1))-.4*R;
%         zshift = .8*R*(rand(obj.lesionFeatures,1))-.4*R;
%         for k = 1:obj.lesionFeatures
%            % create lobule
%            [x,y,z] = sphere(150);
%            x = 0.8*R*x(:) + xshift(k);
%            y = 0.8*R*y(:) + yshift(k);
%            z = 0.8*R*z(:) + zshift(k);
%            x1 = [x1; x];
%            y1 = [y1; y];
%            z1 = [z1; z];
%         end
%         lesionpoints = struct('x',x1,'y',y1,'z',z1);
    elseif type == 3    %irregular
        % add 15 small lobules
        numB = obj.lesionFeatures;
        theta = 2*pi*rand(1,numB);
        phi = asin(2*rand(1,numB)-1);
        [xadd,yadd,zadd] = sph2cart(theta,phi,.95*R*ones(1,numB));
        for k = 1:numB
           [x,y,z] = sphere(150);
           x = obj.lesionFeatureSize*x(:) + xadd(k);
           y = obj.lesionFeatureSize*y(:) + yadd(k);
           z = obj.lesionFeatureSize*z(:) + zadd(k);
           x1 = [x1; x];
           y1 = [y1; y];
           z1 = [z1; z];
        end
        lesionpoints = struct('x',x1,'y',y1,'z',z1);
    elseif type == 4
        % generate spicule points, if desired
        num = obj.lesionFeatures;
        spts = GenerateSpicule([x1,y1,z1],num,R);
        % combine data together
        lesionpoints.x = [x1;spts.x];
        lesionpoints.y = [y1;spts.y];
        lesionpoints.z = [z1;spts.z];
    end
end

%%
function spiculepoints = GenerateSpicule(lpts,num,R)
% make 3D spicule model

spiculepoints = struct('x',[],'y',[],'z',[]);

% base parameters %
% number of surrounding points
pn = 30;
% growth parameters %
% directional scatter
scatdev = .15;
% segment length
lavg = R/50;
% thickness start
t0 = R/3;
%thickness decrease
tdecavg = t0/100;
%thickness scatter
tscat = .75*t0;

for scount = 1:num
% starting point
index = randi(numel(lpts(:,1)));

p0 = .8*[lpts(index,1),lpts(index,2),lpts(index,3)];

%generate points
%preallocate
xps = zeros(500,1);
yps = zeros(500,1);
zps = zeros(500,1);
xts = zeros(10000,1);
yts = zeros(10000,1);
zts = zeros(10000,1);
%initialize growth
xps(1) = p0(1);
yps(1) = p0(2);
zps(1) = p0(3);
idx = 1;
tidx = 1;
t = zeros(1000);
t(1) = t0+tscat*rand(1);
%initialize direction
r = sqrt(sum(p0.^2));
ele = pi/2-acos(p0(3)/r);
azi = atan2(p0(2),p0(1));
th = linspace(0,2*pi,pn);
while t(idx) > 0
   %pick length of segment
   lcur = 2*lavg*rand(1);
   %calculate position of next point
   xps(idx+1) = xps(idx) + lcur*cos(azi);
   yps(idx+1) = yps(idx) + lcur*sin(azi);
   [xdif,ydif,zdif] = sph2cart(azi,ele,lcur);
   xps(idx+1) = xps(idx) + xdif;
   yps(idx+1) = yps(idx) + ydif;
   zps(idx+1) = zps(idx) + zdif;
   
   % calculate surrounding points
   %normal vectors
   n1 = cross([xps(idx),yps(idx),zps(idx)],[-1,-1,-1]);
   n1 = n1./norm(n1);
   n2 = cross(n1,[xps(idx),yps(idx),zps(idx)]);
   n2 = n2./norm(n2);
   for p = 1:pn
       xts(tidx) = xps(idx) + t(idx)*cos(th(p))*n1(1)...
                            + t(idx)*sin(th(p))*n2(1);
       yts(tidx) = yps(idx) + t(idx)*cos(th(p))*n1(2)...
                            + t(idx)*sin(th(p))*n2(2);
       zts(tidx) = zps(idx) + t(idx)*cos(th(p))*n1(3)...
                            + t(idx)*sin(th(p))*n2(3);
       tidx = tidx+1;
   end
   
   %update direction
   ele = ele + scatdev*(rand(1)-.5);
   azi = azi + scatdev*(rand(1)-.5);
   %update thickness
   t(idx+1) = t(idx) - exprnd(tdecavg);
   idx = idx + 1;
end

%truncate
xts = xts(1:tidx-1);
yts = yts(1:tidx-1);
zts = zts(1:tidx-1);



%add points to outvar
spiculepoints.x = [spiculepoints.x;xts];
spiculepoints.y = [spiculepoints.y;yts];
spiculepoints.z = [spiculepoints.z;zts];

end
end

%%
function LesVol = ConvertToVolume(data)
% convert to volume data

% shift and scale data
xmin = min(data.x);
ymin = min(data.y);
zmin = min(data.z);

x = data.x - xmin+.0001;
y = data.y - ymin+.0001;
z = data.z - zmin+.0001;

xr = ceil(max(x));
yr = ceil(max(y));
zr = ceil(max(z));

% create voxel size
dz = max([xr yr zr])/data.res;
dy = dz;
dx = dz;
xdim = ceil(xr/dx);
ydim = ceil(yr/dy);
zdim = ceil(zr/dz);
% create empty volume
LesVol = zeros([xdim,ydim,zdim]);

% check for intersections layer by layer
zstart = 0;
for laynum = 1:zdim
   inds = find(z>zstart&z<=zstart+dz);
   xvals = ceil(x(inds)./dx);
   yvals = ceil(y(inds)./dy);
   layinds = laynum*ones(size(xvals));
   VolInds = sub2ind(size(LesVol),xvals,yvals,layinds);
   LesVol(VolInds) = 1;
   LesVol(:,:,laynum) = bwmorph(LesVol(:,:,laynum),'bridge');
   LesVol(:,:,laynum) = imfill(LesVol(:,:,laynum),'holes');
   zstart = zstart+dz;
end

end
%%
function LesionPlot(data)

figure(9);
p = patch(isosurface(data.lesvol));
set(p,'FaceColor','red');
set(p,'EdgeColor','none');
daspect([1,1,1]);
axis tight
view(-38,30)
camlight
lighting gouraud
title('IsoSurface Representation')
end
