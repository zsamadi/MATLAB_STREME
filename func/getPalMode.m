function [palModel, edOptimum]=getPalMode(PWMB, PWM0)
W=size(PWMB, 2);

bestModel=PWMB;

rvModel=bestModel(:, end:-1:1);

ModelEx=repmat(PWM0,1, 3*W-2);

bestModelEx=ModelEx;
bestModelEx(:, W:2*W-1)=bestModel+PWM0;

edOptimum=norm(bestModelEx);

for shift=-2:2
    rcModelEx=ModelEx;

    rcModelEx(:, W+shift:2*W+shift-1)=rvModel+PWM0;




    eDist=sqrt(norm(bestModelEx-rcModelEx));

    if eDist<edOptimum
        edOptimum=eDist;
        rcModelExOptimum=rcModelEx;
        shiftOptimum=shift;
    end
end

palModel=(bestModelEx+rcModelExOptimum)-2*PWM0;
palModel(:, W+shiftOptimum:2*W+shiftOptimum-1)=palModel(:, W+shiftOptimum:2*W+shiftOptimum-1)/2;

palModel=palModel(:, W+floor(shiftOptimum/2):2*W+floor(shiftOptimum/2)-1);