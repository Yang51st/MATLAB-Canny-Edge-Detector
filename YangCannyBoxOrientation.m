%Loading in the image and getting it grayscaled.
img1=imread("BoxQ.jpeg");
img1=double(img1)/255; %Normalizing the image so that values are between 0 and 1.
gscale=[0.2989,0.5870,0.1140]; %Values for grayscaling an RGB image.
chan1=img1(:,:,1); %Calculating values from the channels of the image.
chan2=img1(:,:,2);
chan3=img1(:,:,3);
img1=gscale(1).*chan1+gscale(2).*chan2+gscale(3).*chan3; %Combining the channels to make the image grayscaled.
rowL=length(img1(:,1)); %Variable for the number of pixel rows in the image.
colL=length(img1(1,:)); %Variable for the number of pixel columns in the image.


%Creating Gaussian filter
filtSize=7.0; %Determines the size of the Gaussian filter.
sigval=1; %Sigma value for the Gaussian filter.
pos=-floor(filtSize/2):floor(filtSize/2); %Creating a grid for the 2D Gauss filter values.
[gfiltx,gfilty]=meshgrid(pos,pos);
gfilt=1./(2.*pi().*sigval.^2).*exp(-1.*(gfiltx.^2.+gfilty.^2)./(2.*sigval.^2)); %Creating 2D Gauss filter.
gfilt=gfilt./sum(gfilt(:)); %Normalizing the Gauss filter.
bimg1=img1;


%Applying horizontal and vertical gradients
gfilth=-1.*gfiltx./sigval.^2.*gfilt; %Creating horizontal derivative-of-Gaussian filter.
gfiltv=-1.*gfilty./sigval.^2.*gfilt; %Creating vertical derivative-of-Gaussian filter.
gfilth=gfilth/sum(sum(abs(gfilth))); %Normalizing.
gfiltv=gfiltv/sum(sum(abs(gfiltv))); %Normalizing.
himg1=imfilter(bimg1,gfilth); %Image with horizontal filter applied.
vimg1=imfilter(bimg1,gfiltv); %Image with vertical filter applied.


%Combining horizontal and vertical gradients to get magnitudes and angles
mimg1=sqrt(himg1.^2.+vimg1.^2); %Image with magnitude of gradients.
aimg1=atand(vimg1./himg1); %Image with angle of gradients.


%Supress the edges that are not maximums along a direction
frameSize=3; %Variable for size of non-max suppression frame size.
vertM=zeros(frameSize,frameSize); %Vertical suppression matrix.
horM=zeros(frameSize,frameSize); %Horizontal suppresion matrix.
diag1M=zeros(frameSize,frameSize); %First diagonal suppression matrix.
diag2M=zeros(frameSize,frameSize); %Second diagonal suppression matrix.
for row=1:frameSize %For-loop for populating the suppression matrices.
    for col=1:frameSize
        if (row==col)
            diag1M(row,col)=1;
        end
        if (col==ceil(frameSize/2))
            vertM(row,col)=1;
        end
        if (row==ceil(frameSize/2))
            horM(row,col)=1;
        end
        if ((frameSize-row+1)==col)
            diag2M(row,col)=1;
        end
    end
end
simg1=zeros(size(mimg1)); %This image, at first black, will be populated with edges that survive non-max suppression.
mimg1=padarray(mimg1,[floor(frameSize/2),floor(frameSize/2)]); %Padding the magnitude of gradients image so that the suppression matrices can be convolved.
for row=ceil(frameSize/2):(rowL+floor(frameSize/2)) %For-loops will go through each pixel, check its gradient angle, then choose the appropriate suppresion matrix.
    for col=ceil(frameSize/2):(colL+floor(frameSize/2))
        checkM=mimg1(row-floor(frameSize/2):row+floor(frameSize/2),col-floor(frameSize/2):col+floor(frameSize/2));
        if (aimg1(row-floor(frameSize/2),col-floor(frameSize/2))>=67.5 && aimg1(row-floor(frameSize/2),col-floor(frameSize/2))<=90)
            dirM=vertM;
        end
        if (aimg1(row-floor(frameSize/2),col-floor(frameSize/2))>=-90 && aimg1(row-floor(frameSize/2),col-floor(frameSize/2))<=-67.5)
            dirM=vertM;
        end
        if (aimg1(row-floor(frameSize/2),col-floor(frameSize/2))>=-22.5 && aimg1(row-floor(frameSize/2),col-floor(frameSize/2))<=22.5)
            dirM=horM;
        end
        if (aimg1(row-floor(frameSize/2),col-floor(frameSize/2))>22.5 && aimg1(row-floor(frameSize/2),col-floor(frameSize/2))<67.5)
            dirM=diag1M;
        end
        if (aimg1(row-floor(frameSize/2),col-floor(frameSize/2))>-67.5 && aimg1(row-floor(frameSize/2),col-floor(frameSize/2))<22.5)
            dirM=diag2M;
        end
        tempM=checkM.*dirM; %Applying the suppression matrix.
        if (max(tempM(:))==tempM(ceil(frameSize/2),ceil(frameSize/2))) %If pixel passes non-max suppression, its value is moved into the simg1 image.
            simg1(row-floor(frameSize/2),col-floor(frameSize/2))=mimg1(row,col);
        end
    end
end


%Applying threshold
higTh=0.05; %Any pixel value lower than this threshold is put to 0, and greater or equal to it is raised to 1.
fimg1=simg1;
for row=1:rowL
    for col=1:colL
        if (fimg1(row,col)>=higTh)
            fimg1(row,col)=1;
        else
            fimg1(row,col)=0;
        end
    end
end %Now the final image, fimg1, contains only pixels with values for either edge or no edge.
imshow(fimg1);


%Using Canny edge detector output to determine orientation of box in image.
edgePixels=0; %Counts the number of edge pixels.
rightPixels=0; %Counts the number of edge pixels that are horizontal or vertical.
tole=10; %Tolerance for the variation in degrees.
orion=0.6; %Threshold percentage at and past which the box is considered properly aligned.
for row=ceil(frameSize/2):(rowL-floor(frameSize/2)) %For-loops that skip the pixels at the image boundary.
    for col=ceil(frameSize/2):(colL-floor(frameSize/2))
        if fimg1(row,col)
            edgePixels=edgePixels+1; %Conditional statements below check if edge is horizontal or vertical.
            if (aimg1(row,col)>=(90-tole) && aimg1(row,col)<=(90+tole))
                rightPixels=rightPixels+1;
            elseif (aimg1(row,col)>=(-90-tole) && aimg1(row,col)<=(-90+tole))
                rightPixels=rightPixels+1;
            elseif (aimg1(row,col)>=(0-tole) && aimg1(row,col)<=(0+tole))
                rightPixels=rightPixels+1;
            end
        end
    end
end
percy=double(rightPixels/edgePixels); %Percentage of edge pixels that are horizontal or vertical.
disp(percy); %Displaying the percentage.
if percy>=orion %Displaying a message on whether the box in the image is aligned correctly or not.
    disp("Box is oriented correctly.");
else
    disp("Box is not oriented correctly.");
end