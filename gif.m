for j=1:1:865
        im=imread(sprintf('%d.bmp',j));
        [I,map]=rgb2ind(im,20); 
        if j==1
           imwrite(I,map,'new.gif','gif', 'Loopcount',inf,'DelayTime',0.02);%FIRST
       else
           imwrite(I,map,'new.gif','gif','WriteMode','append','DelayTime',0.02);
       end
end
    close all;