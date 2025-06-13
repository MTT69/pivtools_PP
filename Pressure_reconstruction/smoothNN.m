function B = smoothNN(A)
% Returns a near-neighbor smoothing of the 2D matrix A

    B = nanmean(cat(3,A,circshift(A,[-1 -1]),circshift(A,[0 -1]), ...
        circshift(A,[1 -1]), circshift(A,[1 0]),circshift(A,[1 1]), ...
        circshift(A,[0 1]),circshift(A,[-1 1]),circshift(A,[-1 0])),3);

end

