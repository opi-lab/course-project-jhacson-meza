from matplotlib import pyplot as plt
import numpy as np
import warnings
import glob
import cv2
import os

warnings.filterwarnings("ignore",".*GUI is implemented.*")

criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.001)
lk_params = {'criteria': (cv2.TERM_CRITERIA_EPS+cv2.TERM_CRITERIA_COUNT, 100,\
                          0.001), 'maxLevel': 20, 'winSize': (21, 21)}


I = glob.glob('target1/*.jpg')

gps = ['-75.46464622800 10.36930439200 32.33273']



if not os.path.exists('output'):
    os.makedirs('output')




print('\n\n\t\t\t\tDetection of target\n\n')
print('> There are {} images for the target\n'.format(len(I)))

img = cv2.imread(I[0])
im1 = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)

p1 = cv2.goodFeaturesToTrack(im1,100000,0.01,10)
p1 = cv2.cornerSubPix(im1,p1,(5,5),(-1,-1),criteria)

plt.figure('Image {}'.format(os.path.basename(I[0])))
plt.axis('off'), plt.imshow(img[...,::-1])
plt.plot(p1[...,0],p1[...,1],'r.',markersize=2), plt.show()
poi = plt.ginput(3,mouse_add=3,mouse_pop=1,mouse_stop=2,timeout=0)
plt.close()

ind1 = (p1[...,0]>poi[0][0]) & (p1[...,0]<poi[1][0]) & \
(p1[...,1]>poi[0][1]) & (p1[...,1]<poi[2][1])

coor = p1[ind1].ravel()

print('Image {} processed: detection 1/{} completed'\
      .format(os.path.basename(I[0]),len(I)))
    
    
for j, image in enumerate(I[1:],2):
    image_name = os.path.basename(image)

    imc = cv2.imread(image)
    im2 = cv2.cvtColor(imc,cv2.COLOR_BGR2GRAY)

    p2_track,*_ = cv2.calcOpticalFlowPyrLK(im1,im2,p1,None,**lk_params)

    p2 = cv2.goodFeaturesToTrack(im2,100000,0.01,10)
    p2 = cv2.cornerSubPix(im2,p2,(5,5),(-1,-1),criteria)

    d = np.linalg.norm(p2-p2_track[ind1],None,2)
    ind2 = np.nonzero(d==d.min())[0][0]

    coor = p2[ind2].ravel()
    
    if d.min() > 1.8:
        print('\x1b[1;31;400m'+'Target lost in the image {}'.format(\
                    image_name)+'\x1b[0m')

        plt.figure('Imagen {}'.format(image_name))
        plt.axis('off'), plt.imshow(imc[...,::-1])
        plt.plot(p2[...,0],p2[...,1],'r.',markersize=2), plt.show()
        poi = plt.ginput(3,mouse_add=3,mouse_pop=1,mouse_stop=2,timeout=0)
        plt.close()

        ind2 = (p2[...,0]>poi[0][0]) & (p2[...,0]<poi[1][0]) & \
        (p2[...,1]>poi[0][1]) & (p2[...,1]<poi[2][1])

        coor = p2[ind2].ravel()


    imwr = imc.copy()
    cv2.circle(imwr, (coor[0],coor[1]), 10, (0,0,255), -1)
    cv2.imwrite('output/{}'.format(image_name),imwr)


    print('Image {} processed: detection {}/{} completed'\
          .format(image_name,j,len(I)))


    im1 = im2.copy()
    p1 = p2.copy()
    ind1 = ind2
