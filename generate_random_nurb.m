function nurb = generate_random_nurb()

    nurb.coords = zeros(3,3,3);
    nurb.weights = zeros(1,3,3);
    nurb.coords(1,:,:) = reshape([0.00000  1.00000  2.00000  0.00000  1.00000  2.00000  0.00000  1.00000  2.00000], 3, 3);
    nurb.coords(2,:,:) = reshape([0.00000  0.00000  0.00000  1.00000  1.00000  1.00000  2.00000  2.00000  2.00000], 3, 3);
    nurb.weights(1,:,:) = reshape([1.0000   1.00000   1.00000   1.00000   1.0000   1.00000   1.00000   1.00000   1.00000], 3, 3);
    nurb.coords(1:2,:,:) = nurb.coords(1:2,:,:) - 0.5 + rand(2,3,3);
    nurb.weights = rand(1,3,3) * 2;
    nurb.coords(1,:,:) = nurb.coords(1,:,:) .* nurb.weights;
    nurb.coords(2,:,:) = nurb.coords(2,:,:) .* nurb.weights;
    nurb.number = [ 3 3 ];
    nurb.order = [ 3 3 ];
    nurb.knots{1} = [0 0 0 1 1 1];
    nurb.knots{2} = [0 0 0 1 1 1];

end