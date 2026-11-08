function iter=drawall_2d(NUM,X)
global omega N
iter=size(X,2);
ib=1:iter; 
path=moviein(length(ib));p=1;
for ii=ib
    drawnlc_2d(X(end-2*N+1:end,ii));
%     drawnlc_beta(X(end-2*N+1:end,ii));
%     set(gcf,'Position',[200 200 400 428]);
%     set(gca,'Position',[0 0.01 1 0.92]); axis equal;
%     textob = [ '$$E=$$ ' num2str(enenlc_2d(X(end-2*N+1:end,ii))    ,'%.5f,')...
%         '\qquad $$\|\nabla E\|=$$ ' num2str(X(1,ii)','%.2e')];
%     textob = [textob(1:end-4) ' ' textob(end-3:end)];
%     if textob(8)=='-'
%         textob = [textob(1:7) '$$-$$' textob(9:end)];
%     end
%     if textob(end-2)=='+' && textob(end)=='0'
%         textob(end-2)='-';
%     end
%     if textob(end-1)=='0'
%         textob(end-1) = [];
%     end
%     text(0, 1.06,textob,'Interpreter','latex','horiz','center','color',[0 0 0],'fontsize',14);
    set(gcf,'color','w'); 
    path(p)=getframe(gcf);    p=p+1;
end
if isreal(NUM)
    NUM=num2str(NUM);
end
video=VideoWriter(['nlc_2d/movie/V_' num2str(omega,'%.0f') '_' NUM],'MPEG-4');
video.FrameRate=12; open(video); writeVideo(video,path);close(video);
end