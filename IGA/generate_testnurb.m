function nurb = generate_testnurb()

% 3-dimensional coordinates are assumed in most functions.
nurb.coords = zeros(3,2,3);
nurb.weights = zeros(1,2,3);
nurb.coords(1,:,:) = reshape([1.00000   2.00000   0.707106781186548   1.414213562373095   0.00000   0.00000], 2, 3);
nurb.coords(2,:,:) = reshape([0.00000   0.00000   0.707106781186548   1.414213562373095   1.00000   2.00000], 2, 3);
nurb.weights(1,:,:) = reshape([1.00000   1.00000   0.707106781186548   0.707106781186548   1.00000   1.00000], 2, 3);
nurb.number = [ 2 3 ];
nurb.order = [ 2 3 ];
nurb.knots{1} = [0 0 1 1];
nurb.knots{2} = [0 0 0 1 1 1];

end
