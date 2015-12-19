%Script for calling the 3-D Finite Difference function.

clear
clc
close all

Alal = 8.418e-5; %in m^2/s
Alss = 4.2e-6; %in m^2/s
ktc = [247 16];
h = 30; %W/(m^2K) Pretty average for air.

Rcp = [ktc(1)/Alal ktc(2)/Alss];
hspal = h*Alal/ktc(1);
hspss = h*Alss/ktc(2);

Alpha = [Alal Alss];
hsp = [hspal hspss];

Thick = 10;
gs = 5;
gz = 5;
x = 0:gs:350;
y = 0:gs:50;
z = 0:gz:Thick;

[X,Y,Z] = ndgrid(x,y,z);

gs2 = gs*1e-3;
gz2 = gz*1e-3;

slice = {'Bottom','Center','Top'};
timestep = [.005 .005];
Lbl = {'Aluminum','Stainless Steel'};
stop = 5*60;
Tamb = 170; %Ambient T in K.
vid = VideoWriter('HeatEquation.avi');
vid.FrameRate = 5;
open(vid)


for g = 1:2
    
    %Initial T = 22C
    Map = zeros(size(Z)) + 22;
    
    Stab = 2*Alpha(g)*timestep(g)*(2/(gs2^2)+1/(gz2^2));
    fprintf('If %.5f is less than 1 its stable. \n',Stab)
    
    st = 1;
    c = 1;
    k = 0;
    marray = 20;
    c2 = 0;
    
    while st
        
        Map = ThreeDFiniteDiff(Map,Rcp(g),ktc(g),timestep(g),gs2,Tamb,gz2,h);
        c = c+1;
        k= k + timestep(g);
        c2(c) = k;
        
        marray(c)=mean(mean(Map(:,:,1)));
        
        if mod(c,250) == 0
            fig = figure;
            set(fig,'Visible','off')
            for w = 1:length(z)
                subplot(length(z),1,w)
                contour(X(:,:,1),Y(:,:,1),Map(:,:,w),'showtext','on');
                hold on
                title(sprintf('%-s Bar after %.2f seconds; slice: %-s',Lbl{g},k,slice{w}))
                xlabel('x')
                ylabel('y')

            end
            fr = getframe(fig);
            writeVideo(vid,fr);
        end
        
%         if mod(c,10000) == 0
%             figure(g)
%             clf
%             for w = 1:length(z)
%                 subplot(length(z),1,w)
%                 contour(X(:,:,1),Y(:,:,1),Map(:,:,w),'showtext','on');
%                 hold on
%                 title(sprintf('%-s Bar after %.0f seconds; slice: %.0f',Lbl{g},k,w))
%                 xlabel('x')
%                 ylabel('y')
%                 drawnow
%             end
%             
%         end
            
            if k>stop
                st = 0;
            end
            
    end
    if g == 1
        m1 = marray;
        cee1 = c2;
    else
        m2 = marray;
        cee2 = c2;
    end
end


close(vid)
figure
plot(cee1,m1)
hold on
plot(cee2,m2)
grid on
xlabel('Time (s)')
ylabel('Avg T of base in C')
title('Base T of Bar in oven')
legend(Lbl,'Location','NorthWest')