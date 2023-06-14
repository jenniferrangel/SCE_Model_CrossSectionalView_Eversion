%First Run ExtractApiBasNucInfo.m to get apical, basal, nuc and order matrices

%The following lines are added to make sure the basal node, apical node and nucleus center positions and order are imported correctly from the ExtractApiBasNucInfo.m
basal_check = basal(:,:);
apical_check = apical(:,:);
nuc_check = nuc(:,:);
order_check = order(:,:);

%Calculate Relative Nuclear Positioning:
%*****************************************************
%Nuclear positioning is defined as d_b / (d_b + d_a) where 
% d_b and d_a is the distance of the nucleus from the basal and apical surfaces, respectively 

d_b = zeros(size(nuc,1),1);
d_a = zeros(size(nuc,1),1);
nuclear_positioning = zeros(size(nuc,1),1);
nuclear_positioning_CellOrder = zeros(size(nuc,1),2);

if (size(basal,1) ~= size(apical,1))
    disp('Check dimension of your matrices')
else
    for i = 1:size(nuc,1)

        d_b(i) = norm(basal(i,:) - nuc(i,:));
        d_a(i) = norm(apical(i,:) - nuc(i,:));

        nuclear_positioning(i) = d_b(i)/(d_b(i) + d_a(i));

        %To keep track of the nuclear positioning and cell order
        nuclear_positioning_CellOrder(i,1) = order(i);
        nuclear_positioning_CellOrder(i,2) = nuclear_positioning(i);
    end
end

%Now to save the information into an Excel spreadsheet:
writematrix(nuclear_positioning,'nuclearPositioning.xlsx')    %just includes the nuclear positioning
writematrix(nuclear_positioning_CellOrder,'nuclearPositioning_CellOrder.xlsx') %includes the nuclear positiong and cell order all in one matrix
