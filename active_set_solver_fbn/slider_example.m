clc; clear; close all;

fig = uifigure('Position',[100 100 600 600]);

% cg = uigauge(fig,'Position',[100 100 120 120]);
ax = uiaxes(fig,'Position',[50,100,500,500]);
sld = uislider(fig,...
               'Position',[200 50 160 3],...
               'ValueChangingFcn',@(sld,event) sliderMoving(event,ax));
sld.Limits = [0,4];


% Create ValueChangingFcn callback
function sliderMoving(event,ax)
x = event.Value;
U1 = linspace(-5,2-x,100);
U2 = linspace(2-x,5,100);
F1 = 0.5*U1.^2 + 1;
F1s = 0.5*U2.^2 + 1;
F2 = 0.5*U2.^2 + ((x+U2)/2).^2;
F2s = 0.5*U1.^2 + ((x+U1)/2).^2;
cla(ax);
hold(ax,'on');
plot(ax,U1,F1,'-b','linewidth',3);
plot(ax,U2,F1s,':b','linewidth',3);
plot(ax,U2,F2,'-r','linewidth',3);
plot(ax,U1,F2s,':r','linewidth',3);
% plot(ax,-5/6,1.0417,'Marker','o')
plot(ax,[2-x,2-x],[0,10],':k');
axis(ax,'equal');
axis(ax,[-5,5,0,10]);

end