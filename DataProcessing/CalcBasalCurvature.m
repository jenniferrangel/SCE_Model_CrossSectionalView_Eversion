vec_l = zeros(size(basal,1),3);
vec_r = zeros(size(basal,1),3);
vec_z = zeros(size(basal,1),1);
orientation = zeros(size(basal,1),1);
for i = 1:size(basal,1)
    if size(basal,1) == 86
        if (i == 1)
            vec_l(i,:) = (-basal(i,:)+basal(size(basal,1),:));
            vec_r(i,:) = (-basal(i,:)+basal(i+1,:));
            crossproduct = cross(vec_l(i,:),vec_r(i,:));
            if (crossproduct >= 0)
                vec_z(i) = 1;
            else
                vec_z(i) = -1;
            end
            continue;
        end
        if (i == size(basal,1))
            vec_l(i,:) = (-basal(i,:)+basal(i-1,:));
            vec_r(i,:) = (-basal(i,:)+basal(1,:));
            crossproduct = cross(vec_l(i,:),vec_r(i,:));
            if (crossproduct >= 0)
                vec_z(i) = 1;
            else
                vec_z(i) = -1;
            end
            continue;
        end
        vec_l(i,:) = (-basal(i,:)+basal(i-1,:));
        vec_r(i,:) = (-basal(i,:)+basal(i+1,:));
        crossproduct = cross(vec_l(i,:),vec_r(i,:));
        if (crossproduct >= 0)
            vec_z(i) = 1;
        else
            vec_z(i) = -1;
        end
    else
        if i >= size(order,1)
            break;
        end
        if min(order) == 0
            order = order + 1;
        end
        if (i == 1)
%             vec_l(i,:) = (-basal(order(i),:)+basal(order(size(basal,1)),:));
%             vec_r(i,:) = (-basal(order(i),:)+basal(order(i+1),:));
%             crossproduct = cross(vec_l(i,:),vec_r(i,:));
%             if (crossproduct >= 0)
%                 vec_z(i) = 1;
%             else
%                 vec_z(i) = -1;
%             end
            continue;
        end
        if (i == size(basal,1))
            vec_l(i,:) = (-basal(order(i),:)+basal(order(i-1),:));
            vec_r(i,:) = (-basal(order(i),:)+basal(order(1),:));
            crossproduct = cross(vec_l(i,:),vec_r(i,:));
            if (crossproduct >= 0)
                vec_z(i) = 1;
            else
                vec_z(i) = -1;
            end
            continue;
        end
        vec_l(i,:) = (-basal(order(i),:)+basal(order(i-1),:));
        vec_r(i,:) = (-basal(order(i),:)+basal(order(i+1),:));
        crossproduct = cross(vec_l(i,:),vec_r(i,:));
        if (crossproduct >= 0)
            vec_z(i) = 1;
        else
            vec_z(i) = -1;
        end
    end
end

% figure(1); plot(basal(:,1),basal(:,2),'-o');
% axis equal;
length = zeros(size(basal,1),3);
for i = 1:size(basal,1)%2:(size(basal,1)-1)
    if (i == 1)
        length(i,1) = norm(basal(i,:)-basal(size(basal,1),:));
        length(i,2) = norm(basal(i,:)-basal(i+1,:));
        length(i,3) = norm(basal(i+1,:)-basal(size(basal,1),:));
        continue;
    end
    if (i == size(basal,1))
        length(i,1) = norm(basal(i,:)-basal(i-1,:));
        length(i,2) = norm(basal(i,:)-basal(1,:));
        length(i,3) = norm(basal(1,:)-basal(i-1,:));
        continue;
    end
    length(i,1) = norm(basal(i,:)-basal(i-1,:));
    length(i,2) = norm(basal(i,:)-basal(i+1,:));
    length(i,3) = norm(basal(i+1,:)-basal(i-1,:));
end

RofCurv = zeros(size(basal,1),1);
Curv = zeros(size(basal,1),1);
for i = 1:size(basal,1)%2:(size(basal,1)-1)
%     %w = (length(i,3)^2 + length(i,2)^2 - length(i,1)^2)/(2*length(i,3)*length(i,2));
%     w = (length(i,1)^2 + length(i,2)^2 - length(i,3)^2)/(2*length(i,1)*length(i,2));
%     ww = acos(w);
%     www = sin(ww);
% %     RofCurv(i) = length(i,1)/(2*www);
%     RofCurv(i) = length(i,3)/(2*www);
%     Curv(i) = sign(vec_z(i))*1/RofCurv(i);

%% Menger Curvature Approach
s = (length(i,1)+length(i,2)+length(i,3))/2;
A = sqrt(s*(s-length(i,1))*(s-length(i,2))*(s-length(i,3)));
Curv(i) = sign(vec_z(i))*4*A/(length(i,1)*length(i,2)*length(i,3));
end

Height = zeros(size(basal,1),1);
if (size(basal,1) ~= size(apical,1))
else
    for i = 1:size(basal,1)
        Height(i) = norm(apical(i,:) - basal(i,:));
    end
end

% Nucleus_Height = zeros(size(basal,1),1);
% if (size(basal,1) ~= size(apical,1))
% else
%     for i = 1:size(basal,1)
% %         if (i >= size(basal,1)-4)
% %             i
% %             nucleus(i,:)
% %             basal(i,:)
% %         end
%         Nucleus_Height(i) = norm(nucleus(i,:)-basal(i,:));
%     end
% end
% Nucleus_Height_Percentage = Nucleus_Height./Height;
% figure(2)
% plot(2:62, Curv(3:63),'-o')%1./RofCurv(3:63), '-o')
% xlabel('Cell ID','FontSize',15);
% title('K_{contract} = 12\alpha, T = 50,000(A.U.)', 'FontSize', 15)
% ylabel('Local Curvature','FontSize',15)
% if (min(Curv(3:63)) < 0 || max(Curv(3:63))<0)
%     axis([1,63,1.1*min(Curv(3:63)),1.1*max(Curv(3:63))]);
% else
%     axis([1,63,0.9*min(Curv(3:63)),1.1*max(Curv(3:63))]);
% end
% 
% figure(3)
% plot(2:62, Height(3:63), '-o')
% axis([1,63,0.9*min(Height(3:63)),1.1*max(Height(3:63))]);
% title('Cell height (\mum), K_{contract} = 12\alpha, T = 50,000(A.U.)','FontSize',15)
% xlabel('Cell ID','FontSize',15);
% ylabel('Cell Height (\mum)','FontSize',15)
% vec_l = zeros(size(basal,1),2);
% vec_r = zeros(size(basal,1),2);
% for i = 2:(size(basal,1)-1)
%     vec_l(i,1) = basal(i-1,1) - basal(i,1);
%     vec_l(i,2) = basal(i-1,2) - basal(i,2);
%     vec_r(i,1) = basal(i+1,1) - basal(i,1);
%     vec_r(i,2) = basal(i+1,2) - basal(i,2);
% end
% 
% angle = zeros(size(basal,1),1);
% absv = zeros(size(basal,1),1);
% norm_l = zeros(size(basal,1),1);
% norm_r = zeros(size(basal,1),1);
% dotp_div_norms = zeros(size(basal,1),1);
% for j = 2:size(basal,1)-1
%     absv(j) = dot(vec_l(j,:),vec_r(j,:));
%     norm_l(j) = norm(vec_l(j,:));
%     norm_r(j) = norm(vec_r(j,:));
%     dotp_div_norms(j) = absv(j)/(norm_l(j)*norm_r(j));
%     angle(j) = acos(dotp_div_norms(j));
% end

