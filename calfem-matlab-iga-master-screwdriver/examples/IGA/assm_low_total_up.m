

% IN_total=[ in_low;INPUT_NORMAL;in_up];
% OUT_total=[ out_low;OUTPUT_NORMAL;out_up];
% Rank=randi(size(IN_total,1),1,size(IN_total,1));
% IN_total=IN_total(Rank,:);
% OUT_total=OUT_total(Rank);
IN_total_add1=[IN_total_2W; INPUT_];
OUT_total_add1=[OUT_total_2W;OUTPUT_];
Rank=randi(size(IN_total_add1,1),1,size(IN_total_add1,1));
IN_total_add1=IN_total_add1(Rank,:);
OUT_total_add1=OUT_total_add1(Rank);