% dataSize = [dim1 dim2 dim3]
% numOfComp = Number of signal components
% noiseNumOfComp = Number of noise components
% noisePercent = Noise norm as a percentage of signal norm 
% type : 0 -> Noise is generated as the sum of "noiseNumComp" outer 
%             products of vectors with scaled elements from the standard 
%             normal distribution.
%        1 -> Same as 0 but the norm of the noise is approximatelly set by 
%             scaling the gaussian vectors so as not to form the whole 
%             tensor. Useful for big tensor data.
%        2 -> Noise is just a tensor with scaled elements from the standard 
%             normal distribution.
function data = artificial_data_generator(dataSize,numOfComp,noiseNumOfComp,noisePercent,type)

data = [];
R=numOfComp;
A = randn(dataSize(1),R);
B = randn(dataSize(2),R);
C = randn(dataSize(3),R);

noise_r = noiseNumOfComp;
A_n = randn(dataSize(1),noise_r);
B_n = randn(dataSize(2),noise_r);
C_n = randn(dataSize(3),noise_r);

if type ==0
    %------- dense exact - with 100 rank-1 factors of Gaussian noise
    data = ktensor(ones(R,1),A,B,C);
    noise =  ktensor(ones(noise_r,1),A_n,B_n,C_n);
    noise = noise*(1/norm(noise))*(norm(data)*noisePercent/100);
    data = data +noise;
elseif type ==1
    %------- dense fast approximate - with 100 rank-1 factors of Gaussian noise
    data = ktensor(ones(R,1),A,B,C);
    A_n = A_n./sqrt(sum(A_n.^2,1))*(norm(data)*noisePercent/100/sqrt(R));
    B_n = B_n./sqrt(sum(B_n.^2,1));
    C_n = C_n./sqrt(sum(C_n.^2,1));
    data = ktensor(ones(R+noise_r,1),[A A_n],[B B_n],[C C_n]);
elseif type ==2
    %------- dense exact with Gaussian noise
    data = ktensor(ones(R,1),A,B,C);
    noise = randn(dataSize);
    noise = noise/norm(tensor(noise))*norm(data)*noisePercent/100;
    data = tensor(data) + noise;
end
