function [epochs]=remove_edge_epochs(epochs, old_segment_length)

%% [epochs]=remove_edge_epochs(epochs, segment_length);
%   
%   The purpose of this function is to remove "epochs" that lie along
%   segment edges based on the original trials for post-processed data. For
%   example, if a theta trough is found within 50 ms of the edge of a
%   trial, the resegmentation to a trough-locked theta is going to cross
%   the trial border and introduce a discontinuity in the new segment. To
%   ensure this doesn't happen, this function takes 'old_segment_length'
%   and makes sure that no segments are passed within half of the
%   'new_segment_length'. 

%% 

boundaries=1:old_segment_length:max(epochs(:,3));
goodTR=false(size(epochs,1),1);

for i=1:size(epochs,1);
    trialtimes=epochs(i,2):epochs(i,3);
    if any(ismember(trialtimes,boundaries))
        goodTR(i)=0;
    else
        goodTR(i)=1;
    end
end

epochs=epochs(goodTR,:);

end