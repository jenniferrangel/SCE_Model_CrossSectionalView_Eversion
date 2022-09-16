vec_l = zeros(size(apical,1),3);
vec_r = zeros(size(apical,1),3);
vec_z = zeros(size(apical,1),1);
orientation = zeros(size(apical,1),1);
for i = 1:size(apical,1)
    if (i == 1)
        vec_l(i,:) = (-apical(i,:)+apical(size(apical,1),:));
        vec_r(i,:) = (-apical(i,:)+apical(i+1,:));
        crossproduct = cross(vec_l(i,:),vec_r(i,:));
        if (crossproduct >= 0)
            vec_z(i) = 1;
        else
            vec_z(i) = -1;
        end
        continue;
    end
    if (i == size(apical,1))
        vec_l(i,:) = (-apical(i,:)+apical(i-1,:));
        vec_r(i,:) = (-apical(i,:)+apical(1,:));
        crossproduct = cross(vec_l(i,:),vec_r(i,:));
        if (crossproduct >= 0)
            vec_z(i) = 1;
        else
            vec_z(i) = -1;
        end
        continue;
    end
    vec_l(i,:) = (-apical(i,:)+apical(i-1,:));
    vec_r(i,:) = (-apical(i,:)+apical(i+1,:));
    crossproduct = cross(vec_l(i,:),vec_r(i,:));
    if (crossproduct >= 0)
        vec_z(i) = 1;
    else
        vec_z(i) = -1;
    end
end

% figure(1); plot(apical(:,1),apical(:,2),'-o');
% axis equal;
length = zeros(size(apical,1),3);
for i = 1:size(apical,1)%2:(size(apical,1)-1)
    if (i == 1)
        length(i,1) = norm(apical(i,:)-apical(size(apical,1),:));
        length(i,2) = norm(apical(i,:)-apical(i+1,:));
        length(i,3) = norm(apical(i+1,:)-apical(size(apical,1),:));
        continue;
    end
    if (i == size(apical,1))
        length(i,1) = norm(apical(i,:)-apical(i-1,:));
        length(i,2) = norm(apical(i,:)-apical(1,:));
        length(i,3) = norm(apical(1,:)-apical(i-1,:));
        continue;
    end
    length(i,1) = norm(apical(i,:)-apical(i-1,:));
    length(i,2) = norm(apical(i,:)-apical(i+1,:));
    length(i,3) = norm(apical(i+1,:)-apical(i-1,:));
end

RofapiCurv = zeros(size(apical,1),1);
apiCurv = zeros(size(apical,1),1);
for i = 1:size(apical,1)%2:(size(apical,1)-1)
%     %w = (length(i,3)^2 + length(i,2)^2 - length(i,1)^2)/(2*length(i,3)*length(i,2));
%     w = (length(i,1)^2 + length(i,2)^2 - length(i,3)^2)/(2*length(i,1)*length(i,2));
%     ww = acos(w);
%     www = sin(ww);
% %     RofapiCurv(i) = length(i,1)/(2*www);
%     RofapiCurv(i) = length(i,3)/(2*www);
%     apiCurv(i) = sign(vec_z(i))*1/RofapiCurv(i);

%% Menger apiCurvature Approach
s = (length(i,1)+length(i,2)+length(i,3))/2;
A = sqrt(s*(s-length(i,1))*(s-length(i,2))*(s-length(i,3)));
apiCurv(i) = sign(vec_z(i))*4*A/(length(i,1)*length(i,2)*length(i,3));
end

Height = zeros(size(apical,1),1);
if (size(apical,1) ~= size(apical,1))
else
    for i = 1:size(apical,1)
        Height(i) = norm(apical(i,:) - basal(i,:));
    end
end

% Nucleus_Height = zeros(size(apical,1),1);
% if (size(apical,1) ~= size(apical,1))
% else
%     for i = 1:size(apical,1)
% %         if (i >= size(apical,1)-4)
% %             i
% %             nucleus(i,:)
% %             apical(i,:)
% %         end
%         Nucleus_Height(i) = norm(nucleus(i,:)-apical(i,:));
%     end
% end
% Nucleus_Height_Percentage = Nucleus_Height./Height;
% figure(2)
% plot(2:62, apiCurv(3:63),'-o')%1./RofapiCurv(3:63), '-o')
% xlabel('Cell ID','FontSize',15);
% title('K_{contract} = 12\alpha, T = 50,000(A.U.)', 'FontSize', 15)
% ylabel('Local apiCurvature','FontSize',15)
% if (min(apiCurv(3:63)) < 0 || max(apiCurv(3:63))<0)
%     axis([1,63,1.1*min(apiCurv(3:63)),1.1*max(apiCurv(3:63))]);
% else
%     axis([1,63,0.9*min(apiCurv(3:63)),1.1*max(apiCurv(3:63))]);
% end
% 
% figure(3)
% plot(2:62, Height(3:63), '-o')
% axis([1,63,0.9*min(Height(3:63)),1.1*max(Height(3:63))]);
% title('Cell height (\mum), K_{contract} = 12\alpha, T = 50,000(A.U.)','FontSize',15)
% xlabel('Cell ID','FontSize',15);
% ylabel('Cell Height (\mum)','FontSize',15)
% vec_l = zeros(size(apical,1),2);
% vec_r = zeros(size(apical,1),2);
% for i = 2:(size(apical,1)-1)
%     vec_l(i,1) = apical(i-1,1) - apical(i,1);
%     vec_l(i,2) = apical(i-1,2) - apical(i,2);
%     vec_r(i,1) = apical(i+1,1) - apical(i,1);
%     vec_r(i,2) = apical(i+1,2) - apical(i,2);
% end
% 
% angle = zeros(size(apical,1),1);
% absv = zeros(size(apical,1),1);
% norm_l = zeros(size(apical,1),1);
% norm_r = zeros(size(apical,1),1);
% dotp_div_norms = zeros(size(apical,1),1);
% for j = 2:size(apical,1)-1
%     absv(j) = dot(vec_l(j,:),vec_r(j,:));
%     norm_l(j) = norm(vec_l(j,:));
%     norm_r(j) = norm(vec_r(j,:));
%     dotp_div_norms(j) = absv(j)/(norm_l(j)*norm_r(j));
%     angle(j) = acos(dotp_div_norms(j));
% end

