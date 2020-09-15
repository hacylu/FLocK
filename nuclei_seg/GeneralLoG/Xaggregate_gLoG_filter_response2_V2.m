function [Ik] = Xaggregate_gLoG_filter_response2_V2(imageGray, largestSigma,smallestSigma,sigmaStep, thetaStep, kernelSize, alpha)


Ik=zeros(size(imageGray,1),size(imageGray,2),pi/thetaStep+1);

FC=zeros(round(kernelSize+1),round(kernelSize+1));  %% circular multi-scale gLoG kernel
for theta = 0: thetaStep : pi-thetaStep
    F=zeros(round(kernelSize+1),round(kernelSize+1));
    for sigmaX = largestSigma : sigmaStep: smallestSigma;
            for sigmaY = sigmaX : sigmaStep: smallestSigma;
                h = -elipLog([round(kernelSize+1),round(kernelSize+1)], sigmaX, sigmaY, theta);
                h=h*(1+log(sigmaX^alpha))*(1+log(sigmaY^alpha));
                if sigmaX==sigmaY
                    FC=FC+h;
                else
                    F=F+h;                          %% elliptical multi-scale gLoG kernels
                end
            end
    end
    k=theta/thetaStep+2;
    Ik(:,:,k)=imfilter(imageGray,F,'replicate');
end
FC=FC./(pi/thetaStep);
Ik(:,:,1)=imfilter(imageGray,FC,'replicate');

end

%% to show three D surface
% surf(F)
% shading interp
% axis off