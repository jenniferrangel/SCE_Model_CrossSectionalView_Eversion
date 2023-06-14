%First Run ExtractApiBasNucInfo.m to get apical, basal, nuc and order matrices

%The following lines are added to make sure the basal node, apical node and nucleus center positions and order are imported correctly from the ExtractApiBasNucInfo.m
basal_check = basal(:,:);
apical_check = apical(:,:);
nuc_check = nuc(:,:);
order_check = order(:,:);

%Now we calculate the cell height:
Height = zeros(size(basal,1),1);
height_info_matrix = zeros(size(basal,1),3);
Height_CellOrder = zeros(size(basal,1),2);

if (size(basal,1) ~= size(apical,1))
else
    for i = 1:size(basal,1)
        Height(i) = norm(apical(i,:) - basal(i,:));

        height_info_matrix(i,1) = i;
        height_info_matrix(i,2) = order(i);
        height_info_matrix(i,3) = Height(i);

        Height_CellOrder(i,1) = order(i);
        Height_CellOrder(i,2) = Height(i);
    end
end

%Now to save the information into an Excel spreadsheet:
writematrix(Height,'Height.xlsx') %just the cell height
writematrix(Height_CellOrder, 'Height_CellOrder.xlsx') %if we want the cell order and cell height saved
