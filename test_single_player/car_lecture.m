clear; clc; close all;

T = 0.1;
RM = 500;
p = 1000;

b = 1/RM*[T*T/2; T];
A = [1 T; 0 1];

for l = 10:10:100
    C = zeros(2,l);
    C(:,1) = b;
    for t = 2:l
        C(:,t) = A*C(:,t-1);
    end
    control = C'/(C*C')*[p;0];
    norms(l/10) = norm(control);
    
    
    control = flip(control);
    writerObj = VideoWriter(strcat('car', num2str(l)),'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    axis tight
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    hold on;
    title(strcat('l = ', num2str(l)));
    X = zeros(2,l+1);
    for t = 2:l+1
        X(:,t) = A*X(:,t-1) + b*control(t-1);
    end

    plot([-5,p+5],[20,20],'-k');
    plot([-5,p+5],[-20,-20],'-k');
    for t = 1:l+1
        pos = [X(1,t), 0];
        width = 6;
        height = 3;
        xLeft = pos(1)-3;
        yBottom = pos(2)-1.5;
        r = rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'b', 'FaceColor', 'r', 'LineWidth', 4);
        xlim([-5, p+5]);
        ylim([-p/2-5, p/2+5]);
        
        frame = getframe;
        writeVideo(writerObj,frame);
        delete(r);
    end
    close(writerObj);
end

figure; 
plot(10:10:100, norms)
xlabel('Timesteps alloted');
ylabel('Norm of control vector');

X = zeros(2,l+1);
control = flip(control);
for t = 2:101
    X(:,t) = A*X(:,t-1) + b*control(t-1);
end


