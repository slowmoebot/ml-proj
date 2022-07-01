clear
clc

load data_23x23.mat

images = zeros(length(data),size(data(1).x,1),size(data(1).x,2),3);
res_images = zeros(length(data),size(data(1).sol,1),size(data(1).sol,2));

for i = 1:length(data)
    images(i,:,:,1)=((data(i).x ./ data(i).w)+0.5)./3;
    images(i,:,:,2)=((data(i).y ./ data(i).w)+0.5)./3;
    images(i,:,:,3)=data(i).w/2;
    
    imwrite(squeeze(images(i,:,:,:)),sprintf("data/train_23x23/%05d.png",i))

    res_images(i,:,:)=data(i).sol;
    imwrite(squeeze(res_images(i,2:end-1,2:end-1)),sprintf("data/solution_23x23/%05d.png",i))
end