function [tSlice] = getSliceT()
gap = 66; %mean interslice gap in ms
tSignal = (0:32)'*gap;
sliceOrder = reshape([(1:17);[(18:33),NaN]],34,1); %interleaved, odd-starting
sliceOrder = sliceOrder(~isnan(sliceOrder));
sliceOrder = [sliceOrder;sliceOrder];
tSlice = tSignal(sliceOrder);
return

