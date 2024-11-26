from cv2 import imread, cvtColor, COLOR_BGR2RGB
import numpy as np
import matplotlib.pyplot as plt

def proc_img(img_path, save_path, color_target_R, color_target_G, color_target_B, T):
	img = imread(img_path)
	img = cvtColor(img, COLOR_BGR2RGB)

	dT = 24*T/img.shape[1] # How much time one pixel represents
	dY = 1/img.shape[0] # How much intensity (y direction) one pixel represents
	t_series = np.zeros(0)
	y_series = np.zeros(0)

	# Look at every column of pixels and filter out ones within color criteria.
	# Keep the mean as the y value of the series.
	for col in range(img.shape[1]):
		pixels = np.arange(img.shape[0])[(img[:,col,0] >= color_target_R[0]) 
										& (img[:,col,0] <= color_target_R[1]) 
										& (img[:,col,1] >= color_target_G[0]) 
										& (img[:,col,1] <= color_target_G[1]) 
										& (img[:,col,2] >= color_target_B[0]) 
										& (img[:,col,2] <= color_target_B[1])]
		if len(pixels) > 0:
			t_series = np.append(t_series, col*dT)
			y_series = np.append(y_series, 1-pixels.mean()*dY)

	np.savetxt(save_path, np.column_stack((t_series, y_series)))

	plt.plot(t_series, y_series)
	plt.xlabel("Time (hours)")
	plt.ylabel("Norm. Intensity")
	plt.xlim(0, 24*T)
	plt.ylim(0, 1)
	plt.show()

img_path = ""
save_path = ""
# Specify the color range of the pixels to select
color_target_R = [150, 255]
color_target_G = [0, 100]
color_target_B = [0, 100]
T = 7 # Number of days represented in the image
proc_img(img_path, save_path, color_target_R, color_target_G, color_target_B, T)

