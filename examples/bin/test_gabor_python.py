import cv2
import numpy as np
import matplotlib.pyplot as plt
from compute_FOV_mask import *
from computeAUCROC import *
from gabor import *
import sys

print('Testing gabor filter')

img_cat = np.zeros([0, 300, 300])
img_gt_cat = np.zeros([0, 300, 300], dtype=np.int8)
img_mk_cat = np.zeros([0, 300, 300], dtype=np.int8)

for i in range(1, 2):
    print('D:/test_data/Angios_134/' + str(i) + '.png')
    img = cv2.imread('D:/test_data/Angios_134/' + str(i) + '.png', 0) / 255.0
    print(img[0, 0])
    img_gt = cv2.imread('D:/test_data/Angios_134/' + str(i) + '_gt.png', 0) / 127.5 - 1.0
    img_gt = img_gt.astype(np.int8)
    #img = fill_blackarea(img, mask, 50)
    #print(img[0, 0])
    img = 1.0 - img
    mask = np.ones(img.shape, dtype=np.int8)
    
    img_cat = np.concatenate((img_cat, img[np.newaxis, :, :]), axis=0)
    img_gt_cat = np.concatenate((img_gt_cat, img_gt[np.newaxis, :, :]), axis=0)
    img_mk_cat = np.concatenate((img_mk_cat, mask[np.newaxis, :, :]), axis=0)


gab_resp = gaborFilter(img_cat, np.array([15]), 2.65, 180, img_mk_cat)

print(img_cat[0, 0, 0])
plt.subplot(1, 3, 1)
plt.imshow(img_mk_cat[0, :, :], 'gray')
plt.subplot(1, 3, 2)
plt.imshow(img_cat[0, :, :], 'gray')
plt.subplot(1, 3, 3)
plt.imshow(gab_resp, 'gray')
plt.show()

print(aucROC(gab_resp, img_gt_cat))
print(gab_resp.shape)

sys.exit()

plt.subplot(2, 2, 1)
plt.imshow(gab_resp[0, :, :], 'gray')
plt.subplot(2, 2, 2)
plt.imshow(gab_resp[1, :, :], 'gray')
plt.subplot(2, 2, 3)
plt.imshow(gab_resp[2, :, :], 'gray')
plt.subplot(2, 2, 4)
plt.imshow(gab_resp[3, :, :], 'gray')
plt.show()
