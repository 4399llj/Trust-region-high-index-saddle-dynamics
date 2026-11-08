function [] = drawnlc_2d(u,NAME)
global I omega x y s
[q1, q2] = vec2q_2d(u);S = 2 * sqrt(q1.^2 + q2.^2);
v1 = sqrt(0.5 + q1./S);v2 = sqrt(0.5 - q1./S) .* (1-2*(q2<0));
set(gcf,'Position',[200 200 400 400]);subplot(1,1,1);
fig = pcolor(x, y, S'/s);caxis([0 1]);shading interp;hold on;
xlim([-1 1]); ylim([-1 1]);
set(fig, 'EdgeAlpha',0);set(fig, 'LineStyle','none');
set(gca,'Position',[0 0 1 1]);
axis off; xlim([-1 1]);  ylim([-1 1]); axis equal;
colormap([[linspace(254,238,5) linspace(233,209,5) linspace(204,164,5) linspace(158,130,5) linspace(123, 36  ,30) linspace(  33, 8.5,50)]' ...
          [linspace(250,234,5) linspace(229,209,5) linspace(204,164,5) linspace(158,138,5) linspace(133, 60.5,30) linspace(  58,33.5,50)]' ...
          [linspace(  1, 21,5) linspace( 34, 86,5) linspace(100,164,5) linspace(158,158,5) linspace(158.5,173,30) linspace(173.5,198,50)]']/256);
% colormap([[linspace(254,238,5) linspace(233,209,5) linspace(204,164,5) linspace(158,130,5) linspace(123, 45  ,40) linspace(44,35,10)   linspace(34.5,  10,50)]' ...
%           [linspace(250,234,5) linspace(229,209,5) linspace(204,164,5) linspace(158,138,5) linspace(133, 59.5,50)                      linspace(  58,33.5,50)]' ...
%           [linspace(  1, 21,5) linspace( 34, 86,5) linspace(100,164,5) linspace(158,158,5) linspace(158.5,173,30) linspace(173.5,183,20) linspace(183.5,208,50)]']/256);
fig = gcf;fig.PaperPositionMode = 'auto';fig_pos = fig.PaperPosition;fig.PaperSize = [fig_pos(3) fig_pos(4)];

rng(0); str=(S(2:I+1,2:I+1)).^0.25; len=S(2:I+1,2:I+1)+s/10;
rnd=double((randn(I)./str)>1.63);rnd(rnd==0)=NaN;
quiver(x(2:I+1), y(2:I+1),  0.25 * v1(2:I+1,2:I+1)'.* rnd'.*len',...
    0.25 * v2(2:I+1,2:I+1)'.* rnd'.*len', 3,  'color',[0.99 0.99 0.99],...
    'ShowArrowHead', 'off','LineWidth',1);
quiver(x(2:I+1), y(2:I+1), -0.25 * v1(2:I+1,2:I+1)'.* rnd'.*len',...
    -0.25 * v2(2:I+1,2:I+1)'.* rnd'.*len', 3,  'color',[0.99 0.99 0.99],...
    'ShowArrowHead', 'off','LineWidth',1);
if nargin==2
%     print(gcf,sprintf('./nlc_2d/solutions/N_%.0f_%s',omega,NAME),'-depsc','-r288','-opengl');
    print(gcf,sprintf('./nlc_2d/solutions/N_%.0f_%s',omega,NAME),'-dpng','-r288','-opengl');
%     print(gcf,sprintf('./nlc_2d/solutions/N_%.0f_%s',omega,NAME),'-dpdf');
end
hold off;

end