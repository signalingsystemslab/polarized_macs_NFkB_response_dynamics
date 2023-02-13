function [responder_index, responders_fraction, off_times] = get_activity_metrics(smoothed_by_sigma, Wliml, Wlimu, OnThresh, blockLengthThresh)
%Determine which cells are responders and compute off-times
% adapted from Luecke S
NFkBBySigma = smoothed_by_sigma(:,Wliml:Wlimu);
block_length_readout = zeros(size(NFkBBySigma,1),100);
block_length_readout(:,:) = 111;
responder_index = nan(size(smoothed_by_sigma, 1), 1);
for jj = 1:size(NFkBBySigma,1)
        thresh_idx = diff(NFkBBySigma(jj,:)>OnThresh,1,2);
        thresh_start = find(thresh_idx == 1);
        thresh_stop = find(thresh_idx == -1);
        if any(NFkBBySigma(jj,:)>OnThresh,2)
            if isempty(thresh_start) && isempty(thresh_stop)
                block_length = Wlimu-Wliml;
            elseif isempty(thresh_start)
                block_length = thresh_stop;
            elseif isempty(thresh_stop)
                block_length = numel(thresh_idx)+1 - thresh_start;
            elseif ~isempty(thresh_start) && ~isempty(thresh_stop)
                if (thresh_start(1)<thresh_stop(1)) && (thresh_stop(end)>thresh_start(end))
                       block_length = thresh_stop - thresh_start;
                elseif thresh_start(1)<thresh_stop(1) && thresh_start(end)>thresh_stop(end)
                       block_length = [thresh_stop - thresh_start(1:end-1),numel(thresh_idx)+1 - thresh_start(end)];
                elseif thresh_stop(1)<thresh_start(1) && thresh_stop(end)>thresh_start(end)
                       block_length = [thresh_stop(1),thresh_stop(2:end)-thresh_start];
                elseif thresh_stop(1)<thresh_start(1) && thresh_start(end)>thresh_stop(end)
                       block_length = [thresh_stop(1),thresh_stop(2:end)-thresh_start(1:end-1), numel(thresh_idx)+1 - thresh_start(end)];
                end                    
            end        
        else
                block_length = 0;
        end 
        responder_index(jj,1) = any(block_length>=blockLengthThresh, 2);
        block_length_readout(jj, 1:numel(block_length)) = block_length;
end
responders_fraction = nnz(responder_index)/numel(responder_index);

%determine off times
on_array = zeros(size(smoothed_by_sigma(:,Wliml:end)));
off_times = zeros(size(smoothed_by_sigma,1),1);
for ii = 1:size(smoothed_by_sigma,1)
    n = 0;
    for jj = Wliml:size(smoothed_by_sigma,2)
        if smoothed_by_sigma(ii,jj)> OnThresh
            n = n+1;
        else
        n = 0;
        end
        on_array(ii,(jj-Wliml+1)) = n;
    end
    if find(on_array(ii,:)==1, 1)> Wlimu %ignore cells activating for first time after expected activity window
       off_times(ii) = 0;
    else
        if ~isempty(find(on_array(ii,:)>= blockLengthThresh, 1, 'last'))
           off_times(ii) = find(on_array(ii,:)>=blockLengthThresh, 1, 'last');
        else
           off_times(ii) = 0;
        end
    end
end
end

