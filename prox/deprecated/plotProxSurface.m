% Create surface from proximal points

dProxPoints = testPoints-proxPoints;

%Generate the surface
DT = delaunayTriangulation(dProxPoints');
% TR = triangulation(DT.ConnectivityList,proxPoints')
% 
% %Create and apply the colormap
% c = zeros(64,3);
% for i = 1:64
%     c(i,1) = .5;%(i+32)/100;
%     c(i,2) = 0;%i/100;
%     c(i,3) = 0;%i/100;
% end
% colormap(c);
%Plot the surface
figure;
tetramesh(DT,'EdgeColor','none','FaceColor','b','FaceAlpha',0.1)
grid on
box on
xlabel('eta')
ylabel('y_1')
zlabel('y_2')
xlim([-const 0])
ylim([-const/4 const/4])
zlim([-const/4 const/4])
set(gca,'FontSize',30)

hold on
for i=1:n
    caseInd = caseMat(1,i);
    plot3([testPoints(1,i),dProxPoints(1,i)],[testPoints(2,i),dProxPoints(2,i)],[testPoints(3,i),dProxPoints(3,i)],...
        [colorVec{caseInd},'-'],'LineWidth',2)
end     
