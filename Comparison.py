from __future__ import division
import os
import shutil
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import SimpleITK as sitk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfilename
from matplotlib import colors
from Tkinter import Tk

###
startslice = 290
endslice = 300


height = -1
width = -1

colormap = []

tp = 0 #true positive
tn = 0 #true negative
fp = 0 #false positive
fn = 0 #false negative
tani = 0


def compare(knnmrf, groundtruth):
	global colormap, tp, tn, fp, fn

	for i in range(height):
		for j in range(width):
			if knnmrf[i][j] != 0 and groundtruth[i][j] == 1:
				colormap[i][j] = 3 #white for common area
				tp += 1
			elif groundtruth[i][j] == 1:
				colormap[i][j] = 2 #green for ground truth
				fn += 1
			elif knnmrf[i][j] != 0:
				colormap[i][j] = 1 #cyan for knn mrf
				fp += 1
			else:
				colormap[i][j] = 0
				tn += 1


def openMRI(title):
	Tk().withdraw()
	filename = askopenfilename(title = title) 
	image = sitk.ReadImage(str(filename))

	return filename, sitk.GetArrayFromImage(image)


def writeStats(folder):
	textfile = open(folder + "\\Stats.txt", "w")
	textfile.write("True positive: %s\n" % tp)
	textfile.write("True negative: %s\n" % tn)
	textfile.write("False positive: %s\n" % fp)
	textfile.write("False negative: %s\n" % fn)
	textfile.write("Tanimoto Coefficient: " + str(tani) + "\n")
	textfile.close()


def computeTani():
	return tp/(tp+fp+fn)


if __name__ == "__main__" :
	print "> Load KNN-MRF segmented image."
	directory1, knnmrf = openMRI('Load KNN-MRF segmented image')
	print "\tLoaded ", directory1

	print "> Load corresponding ground truth."
	directory2, groundtruth = openMRI('Load corresponding ground truth')
	groundtruth = groundtruth[startslice:endslice+1]
	print "\tLoaded ", directory2

	height = len(groundtruth[0])
	width = len(groundtruth[0][0])

	folder = directory1.split("/")
	folder = folder[len(folder)-2]
	folder = "SegmentedMRI_Compared\\" + folder
	if os.path.exists(folder): shutil.rmtree(folder)
	try: os.makedirs(folder)
	except OSError as exc: raise

	print "> Start comparison. [ from slice", startslice, "to", endslice, "]"
	for i in range(len(knnmrf)):
		print "\tProcessing slice#", i+startslice, "...\r",
		colormap = [[0 for j in range(width)] for k in range(height)]
		compare(knnmrf[i], groundtruth[i])

		fig, ax = plt.subplots()
		frame = plt.gca()
		frame.axes.get_xaxis().set_visible(False)
		frame.axes.get_yaxis().set_visible(False)

		cmap = colors.ListedColormap(['black', 'deepskyblue', 'gold', 'white'])
		cax = ax.imshow(np.array(colormap), cmap=cmap, vmin=0, vmax=3)
		bar = fig.colorbar(cax, ticks=[0, 1, 2, 3])
		bar.ax.set_yticklabels(['Non-Lesion', 'KNN-MRF', 'Manual', 'Common'])

		filename = folder + "\\Slice" + str(i+startslice) + ".png"
		if os.path.isfile(filename): os.remove(filename)		
		plt.savefig(filename)
		plt.close()

	print "\tProcessing slice#", endslice, "...", 
	print "done"
	
	tani = computeTani()
	writeStats(folder)
	print "> Stats.txt created."
	print ">> Done! Your segmented MR image is compared."

