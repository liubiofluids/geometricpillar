a=1; % lateral size of lattice (normalized)
r0=0.375; % pillar radius for r/a = 15/40
thetmax=acos(2*r0/a);   % cutoff theta
delr=.02;               % stepsize 
epr=0.05;              % distance criteria for residency on pillar surfaces
dthet=delr/(r0+epr);    % angle step for circulation movements corresponding to the linear stepsize 'deltr'
deln=0.1;                 % white noise amplitude of theta
Nmax=400;                 % number of iterations
Tmax=400;               % max path length for 1 trajectory
nmsdstep=200;           % step# for mean-squared displacement (MSD)
msdcum=zeros(Nmax, nmsdstep); % MSD data for Nmax trials
msdedge=delr*exp(linspace(0, log(Tmax/delr), nmsdstep+1)); %edge of the elapsed time/path length for binning
gamma_minus=1.1*0.5;    %range of the angular position of the attractive zone: lower bound
gamma_plus=pi/2-1.1*0.5; %range of the angular position of the attractive zone: upper bound

vdgammaavg=zeros(Nmax, 2); % list of residency path lengths for computing their mean
vp0=linspace(0.45, 1., Nmax); % list of escaping probabilities for plotting

%defining color map
gclmap=jet(64); 
gclmap=gclmap./repmat(sqrt(sum(gclmap'.^2))', 1, 3);
gclmap=gclmap.*repmat(linspace(0.5, 1, length(gclmap))', 1, 3);
gclmap=gclmap+repmat(linspace(0., 0.2, length(gclmap))', 1, 3);
gclmap(gclmap>1)=1;
%end of defining color map

	
nTrajSave=10;
idsel=linspace(0, Nmax, nTrajSave+2); % index of example trajectories to be plotted
idsel=floor(idsel(2:end-1));
dataTraj=nan(length(idsel)*2,  floor(Tmax/delr)+2);
nplt=floor(size(dataTraj, 2)/5);
startarray=zeros(Nmax,2);

for iter=1:Nmax
tim=0; % reset dimensionless time or the path length
intx=0; %initial x
inty=0; %initial y
p0=vp0(iter);
vdgamma=[];
while (intx^2+inty^2< r0^2) % generate a random initial position in the gap among pillars (in the first quadrant)
   intx=0.5*rand;
   inty=0.5*rand;
end

% generate pillar arrays
startarray(iter,1)=intx;
startarray(iter,2)=inty;
xmin=-10;
xmax=10;
ymin=-10;
ymax=+10;
centers=[];
radii=[];
for x=xmin:xmax
    for y=ymin:ymax;
        centers=[centers;[x y]];
        radii=[radii;r0];
    end
end
ind = ceil(rand * size(centers,1));
r=[intx inty]; %position
thet=rand*2*pi; %direction of motion
n=[cos(thet), sin(thet)]; %unit vector of orientation
rcum=r; %cumulative list of r
tcum=tim; %cumulative list of tim
while tim<Tmax
%check pillar collision
p1=[floor(r(1)), floor(r(2))];
p2=[ceil(r(1)), ceil(r(2))];
p3=[floor(r(1)), ceil(r(2))];
p4=[ceil(r(1)), floor(r(2))];
p=[p1; p2; p3; p4];%bottom left, top right, top left, bottom right
coll=0;
for i=1:4
    if (p(i,1)-r(1))^2+(p(i,2)-r(2))^2<=(r0+.9*epr)^2;
        coll=1; % collided
        pill=i; % collided pillar index
    end
end
if coll==1
    thetc=atan2((r(2)-p(pill,2)),(r(1)-p(pill,1)));%find collided angle
    d=[cos(thetc), sin(thetc)];%unit vector to bacterial position from pillar center
    dz=[d,0];
    nz=[n,0];
    c=cross(dz,nz); %cross product to check clockwise or ccw
    dang=0;
    b_escape=rand(1)>p0; % if not escape
    b_continue=mod(thetc, pi/2)<gamma_minus | mod(thetc, pi/2)>gamma_plus; % if is within the attractive zone
    
    if b_escape 
        if ~b_continue
            b_hop=true;
        else
            b_hop=false;
        end
        b_continue=true; % always circulate the pillar
    else
        b_continue=mod(thetc, pi/2)<gamma_minus | mod(thetc, pi/2)>gamma_plus; % if is within the attractive zone
        b_hop=false;
    end
    thetc0=thetc;
    dgamma=0;
    if ~b_continue
	dp=[cos(thetc), sin(thetc)];
	r=p(pill,:)+(r0+epr).*dp;
    end
    while b_continue
    if sign(c(3))>0 
        thetc=thetc+dthet;%go counter clockwise
    else
        thetc=thetc-dthet;%go clockwise
    end

    dp=[cos(thetc), sin(thetc)];
    rp=p(pill,:)+(r0+epr).*dp;% new position
    np=sign(c(3))*[-sin(thetc), +cos(thetc)];%new tangent direction
    tim=tim+(r0+epr)*dthet;
    dang=dang+abs(dthet);
    r=rp;
    n=np;
    rcum=[rcum;r];%position array 
    tcum=[tcum;tim];%time array
    b_continue=mod(thetc, pi/2)<gamma_minus | mod(thetc, pi/2)>gamma_plus; % if is within the attractive zone

    if b_hop
        if b_continue
            b_hop=false; %move out of the repulsive zone
        else
            b_continue=true; %continue moving in the repulsive zone
        end
    else
        if b_continue
            % continue
        else
            b_escape=rand(1)>p0;
            if b_escape
                b_hop=true; %continue moving in the repulsive zone
                b_continue=true;
            else
                b_hop=false;
            end
        end
    end
    dgamma=thetc-thetc0;
    end
    if dgamma~=0
        vdgamma=[vdgamma; dgamma];
    end
else
    rp=r+delr*n;
    np=n+deln*2*(rand(2,1)'-0.5);
    tim=tim+delr;
    r=rp;
    n=np;
    rcum=[rcum;r];
    tcum=[tcum;tim];
end
end

%analyze MSD
deltat=pdist(tcum);
deltar=pdist(rcum);
[Y]=discretize(deltat, msdedge);
for imsd=1:length(msdedge)-1
	indt=Y==imsd;
	msdcum(iter, imsd)=mean(deltar(indt).^2);
end
vdgammaavg(iter, :)=[mean(abs(vdgamma)), std(abs(vdgamma))];
indsave=find(idsel==iter);
if ~isempty(indsave)
	nstpsave=min([size(rcum, 1), size(dataTraj, 2)]);
	dataTraj((indsave-1)*2+1:indsave*2, 1:nstpsave)=rcum(1:nstpsave, :)';
end
end

for ii=1:length(idsel)
vxyctr(ii, :)=floor(nanmean([dataTraj((ii-1)*2+1,1:nplt); dataTraj(2*ii,1:nplt)]'))+floor((rand(1, 2)-.5)*10+.5);
end	


vgcf(3)=figure;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log'); 
hold on;
msdave=mean(msdcum);%take mean over trajectories
msdstd=std(msdcum);%take standdard deviation


%draw example trajectories and pillars
dataTrajx=dataTraj(1:2:end-1, 1:nplt);
dataTrajy=dataTraj(2:2:end, 1:nplt);
vgcf(1)=figure;
delete(findall(gcf, 'Type', 'Line'));

figure(vgcf(1))
hold on;
for i=1:size(dataTrajx, 1)
h=plot(dataTrajx(i, :)-vxyctr(i,1), dataTrajy(i, :)-vxyctr(i,2), '.', 'MarkerSize', 6, 'Color', gclmap(floor((vp0(idsel(i))-min(vp0(:)))/(max(vp0(:))-min(vp0(:)))*(size(gclmap, 1)-1))+1,:));

end

xyrange=[min(dataTrajx(:)), max(dataTrajx(:)), min(dataTrajy(:)), max(dataTrajy(:))];

[xp, yp]=meshgrid(floor(xyrange(1)-9.5:xyrange(2)+10.5), floor(xyrange(3)-9.5:xyrange(4)+10.5));
centers_p=[xp(:), yp(:)];
radii_p=r0*ones(numel(xp), 1);
figure(vgcf(1))
viscircles(centers_p, .95*radii_p,'Color', ones(1,3)*.7);
gpos =[200 200   515   548];
figure(vgcf(1))
set(vgcf(1), 'Position', gpos);
axis equal; 
axis([-1, 1, -1, 1]*9);
set(gca, 'FontSize', 20);
 
 
%plot MSD
vgcf(2)=figure;
figure(vgcf(2));
errorbar(.5*(msdedge(2:end)+msdedge(1:end-1)), msdave, msdstd, 'o');
set(gca, 'XScale', 'log', 'YScale', 'log');


