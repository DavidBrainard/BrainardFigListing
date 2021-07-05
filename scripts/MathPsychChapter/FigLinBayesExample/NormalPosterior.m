function [u,K] = NormalPosterior(data,prioru,priorK,R,noiseK)% This formula derived in Brainard notes, 9/17/92invnoiseK = inv(noiseK);invpriorK = inv(priorK);factor1 = R'*invnoiseK*R;factor2 = inv(invpriorK + factor1);W = factor2*R'*invnoiseK;b = factor2*invpriorK*prioru;u = W*data+b;% This is how I read the same notes in 2021.% The answer comes out the same.M = priorK*R'*inv(R*priorK*R' + noiseK);u1 = M*data + prioru - M*R*prioru; if (max(abs(u-u1)) > 1e-8)    error('Two ways of computing the prior do not agree');end% Covariance matrix of posteriorK = inv(invpriorK + R'*invnoiseK*R);