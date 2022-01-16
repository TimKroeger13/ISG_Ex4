function hsky = skyplot(azim,elev,line_style)
if nargin == 2 , line_style = '*'; end
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;
tc = get(cax,'xcolor');
if ~hold_state
hold on;
zenmax = max(90-elev(:)); zenmax = 15*ceil(zenmax/15);
elmax = 90;
az = 0:pi/50:2*pi;
xunit = sin(az);
yunit = cos(az);
for i=[30 60]
plot(xunit*i,yunit*i,'-','color',tc,'linewidth',1);
end
i=90; plot(xunit*i,yunit*i,'-','color',tc,'linewidth',2);
for i=[15:30:75 105:15:zenmax]
plot(xunit*i,yunit*i,':','color',tc,'linewidth',1);
end
for i=30:30:zenmax
text(0,i,[' ' num2str(90-i)],'verticalalignment','bottom');
end
az = (1:6)*2*pi/12; caz = cos(az); saz = sin(az);
ca = [-caz; caz]; sa = [-saz; saz];
plot(elmax*ca,elmax*sa,'-','color',tc,'linewidth',1);
if zenmax > elmax
plot(zenmax*ca,zenmax*sa,':','color',tc,'linewidth',1);
end
rt = 1.1*elmax;
for i = 1:length(az)
loc1 = int2str(i*30);
if i == length(az)
loc2 = int2str(0);
else
loc2 = int2str(180+i*30);
end
text( rt*saz(i), rt*caz(i),loc1,'horizontalalignment','center');
text(-rt*saz(i),-rt*caz(i),loc2,'horizontalalignment','center');
end
view(0,90);
axis(max(zenmax,elmax)*[-1 1 -1.1 1.1]);
set(cax,'position',[.05 .05 .9 .9])
end
yy = (90-elev).*cos(azim/180*pi);
xx = (90-elev).*sin(azim/180*pi);
q = plot(xx,yy,line_style);
if nargout > 0, hsky = q; end
if ~hold_state, axis('equal'); axis('off'); end
if ~hold_state, set(cax,'NextPlot',next); end
end