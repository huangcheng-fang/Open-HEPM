%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%By Fang Huangcheng @BJTU 
%Last update @2023/9/18
%Email: valy_f@bjtu.edu.cn
function Get_gif(filename,dt,new)
filename=['Output-Files\',filename];
frame = getframe(gcf());
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if new
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',dt);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',dt);
end
end