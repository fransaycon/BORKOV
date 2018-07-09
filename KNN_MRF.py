import gc
import math
import multiprocessing
import os
import sys
import time
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import scipy.misc
import shutil
import SimpleITK as sitk
import tkFileDialog as fd
import Tkinter as tk
from Algorithms import Algorithms
from sklearn.neighbors import KNeighborsClassifier
#from Validation import Validation

mask_const = [[250, 97, 430, 106, 381],
[251, 97, 430, 106, 381],
[252, 97, 430, 106, 381],
[253, 97, 429, 106, 381],
[254, 97, 429, 106, 381],
[255, 97, 429, 106, 381],
[256, 96, 429, 105, 381],
[257, 96, 428, 105, 380],
[258, 96, 428, 105, 380],
[259, 96, 428, 105, 380],
[260, 96, 428, 105, 380]]

#user edit
startslice = 250 #0 indexing
endslice = 260

#params
k = -1
threshold1 = -1
MAPiter = -1
threshold2 = -1

#MRI dimension
slices = -1
height = -1
width = -1

#MRF
x = []
xa = []
xb = []
energy = 0.0


def fileChooser(title):
	root = tk.Tk()
	root.withdraw()
	file_path = fd.askopenfilename(title = title)
	
	return str(file_path)


#MRF CODE - START
def execute(minY, maxY, minX, maxX):
		for i1 in range(minY, maxY+1):
			for i2 in range(minX, maxX+1):
				MAP(i1, i2)


#Maximum a posteriori
def MAP(i1, i2):
	global x, energy
	prob = x[i1][i2]
	i = .01

	ef = energyFunction(prob, i1, i2)
	ef1 = ef
	temp = ef1
	boolean = True

	if prob+i <= 1:
		ef2 = energyFunction(prob+i, i1, i2)
		if ef2 < temp: boolean = False

		while ef2 < temp:
			ef = ef2
			x[i1][i2] = prob+i
			temp = ef2
			i += .01
			if prob+i <= 1: ef2 = energyFunction(prob+i, i1, i2)
			else: break

	if prob-i > 0 and boolean:
		ef2 = energyFunction(prob-i, i1, i2)

		while ef2 < temp:
			ef = ef2
			x[i1][i2] = prob-i
			temp = ef2
			i += .01
			if prob-i > 0: ef2 = energyFunction(prob-i, i1, i2)
			else: break
	
	energy += ef


def energyFunction(prob, i1, i2):
	summation = 0
	neighbors = getNeighbors(i1, i2)

	for j in neighbors:
		summation += affinityFunction(prob, x[j[0]][j[1]])

	if xa is not None: summation += affinityFunction(prob, xa[i1][i2])
	if xb is not None: summation += affinityFunction(prob, xb[i1][i2])

	return summation


def affinityFunction(xi, xj): 
	xtemp = xi
	'''
	bias = (xj-xi)/2
	if xi < threshold and xj < threshold: xtemp = xi+bias
	if (xi < threshold and xj > threshold) or (xi > threshold and xj < threshold):
		xtemp = xi-bias
		if xtemp > 1: xtemp = 1
		if xtemp < 0: xtemp = 0
	'''
	#Bhattacharya coeff
	bc = math.sqrt(xtemp*xj) + math.sqrt((1-xtemp)*(1-xj))
	if bc == 0: return 3
	#Bhattacharya distance
	bd = -math.log(bc)

	return bd


#i1 = row, i2 = column
def getNeighbors(i1, i2):
	index = []

	if i1-1 >= 0:
		index.append([i1-1, i2])
		if i2-1 >= 0: index.append([i1, i2-1])
		if i2+1 < len(x[i1]): index.append([i1, i2+1])

	if i1+1 < len(x): index.append([i1+1, i2])

	return index
#MRF CODE END

#OPTIONS 
#  1 - FLAIR Intensity Only 
#  2 - FLAIR and T1 Intensity Only 
#  3 - x,y,z and FLAIR Intensity 
#  4 - x,y,z, FLAIR and T1 Intensity Only 
def knn(flair,t1,mask,dataset,algorithms):	
	#Initialize needed variables
	start_time = time.clock()
	print "> Start KNN initialization."
	technique = dataset[2] 
	param_1 = dataset[3] 
	param_2 = dataset[4] 
	probabilities = [ [ [0.0 for x in range(width)] for y in range(height)] for z in range(slices) ]
	classifier = KNeighborsClassifier(n_neighbors = k, n_jobs = multiprocessing.cpu_count(),algorithm = 'auto')
	classifier.fit(dataset[0],dataset[1])
	print "\tInitialization is successful "  + str(time.clock() - start_time) + " seconds."

	#Get probabilities using knn          
	start_time = time.clock() 
	print "> Start KNN. [ from slice", startslice, "to", endslice, "]"
	for i in range(len(mask)):		
		print "\tProcessing slice#", mask[i][0], "...\r",
		input = []
		
		for y in range(mask[i][1],mask[i][2]+1):
			for x in range(mask[i][3],mask[i][4]+1):	
					if len(param_1) == 1:	
						if technique == 'N':
							input.append([algorithms.normalize(flair[mask[i][0]][y][x],param_1[0],param_2[0])])
						elif technique == 'S':
							input.append([algorithms.standardize(flair[mask[i][0]][y][x],param_1[0],param_2[0])])
					
					elif len(param_1) == 2: 
						if technique == 'N':
							input.append([algorithms.normalize(flair[mask[i][0]][y][x],param_1[0],param_2[0]),algorithms.normalize(t1[mask[i][0]][y][x],param_1[1],param_2[1])])
						elif technique == 'S':
							input.append([algorithms.standardize(flair[mask[i][0]][y][x],param_1[0],param_2[0]),algorithms.standardize(t1[mask[i][0]][y][x],param_1[1],param_2[1])])
					
					elif len(param_1) == 4:
						if technique == 'N':
							norm_z = algorithms.normalize(mask[i][0],param_1[0],param_2[0])
							norm_y = algorithms.normalize(y,param_1[1],param_2[1])
							norm_x = algorithms.normalize(x,param_1[2],param_2[2])
							norm_FLAIR = algorithms.normalize(flair[mask[i][0]][y][x],param_1[3],param_2[3])
							input.append([norm_z,norm_y,norm_x,norm_FLAIR])
						elif technique == 'S':
							std_z = algorithms.standardize(mask[i][0],param_1[0],param_2[0])
							std_y = algorithms.standardize(y,param_1[1],param_2[1])
							std_x = algorithms.standardize(x,param_1[2],param_2[2])
							std_FLAIR = algorithms.standardize(flair[mask[i][0]][y][x],param_1[3],param_2[3])
							input.append([std_z,std_y,std_x,std_FLAIR])
					
					elif len(param_1) == 5:
						if technique == 'N':
							norm_z = algorithms.normalize(mask[i][0],param_1[0],param_2[0])
							norm_y = algorithms.normalize(y,param_1[1],param_2[1])
							norm_x = algorithms.normalize(x,param_1[2],param_2[2])
							norm_FLAIR = algorithms.normalize(flair[mask[i][0]][y][x],param_1[3],param_2[3])
							norm_T1 = algorithms.normalize(t1[mask[i][0]][y][x],param_1[4],param_2[4])
							input.append([norm_z,norm_y,norm_x,norm_FLAIR,norm_T1])
						elif technique == 'S':
							std_z = algorithms.standardize(mask[i][0],param_1[0],param_2[0])
							std_y = algorithms.standardize(y,param_1[1],param_2[1])
							std_x = algorithms.standardize(x,param_1[2],param_2[2])
							std_FLAIR = algorithms.standardize(flair[mask[i][0]][y][x],param_1[3],param_2[3])
							std_T1 = algorithms.standardize(t1[mask[i][0]][y][x],param_1[4],param_2[4])
							input.append([std_z,std_y,std_x,std_FLAIR,std_T1])	

		output = classifier.predict_proba(input)
		index = 0 
		
		for y in range(mask[i][1],mask[i][2]+1):
			for x in range(mask[i][3],mask[i][4]+1):
				probabilities[mask[i][0]][y][x] = output[index][1]
				index+=1
	
	print "\tProcessing slice#", endslice, "...",
	print "done"
	print "\tKNN is successful "  + str(time.clock() - start_time) + " seconds."
	
	return probabilities


def parse(string):
	substring = string.split("/")
	substring = substring[len(substring)-1]
	substring = substring.split(".")

	return substring[0]


def knn_mrf():
	#mask - Brain Mask 
	#flair - flair array 
	#dataset - dataset array 
	#k - number of neighbors in knn 
	#threshold - for binary segmentation of MRF 
	
	global slices, height, width, endslice

	#Load algorithms module
	algorithms = Algorithms()
	
	#Load Dataset - extract features, labels
	start_time = time.clock()
	print "> Load training data."
	dataset = fileChooser('Load training data')
	features = parse(dataset)
	dataset = algorithms.readTrainingData(dataset)
	print "\tDataset extraction successful. " + str(time.clock() - start_time) + " seconds."
	
	#Load FLAIR MRI
	start_time = time.clock() 
	print "> Load FLAIR input."
	flair = fileChooser('Load FLAIR input')
	mricase = parse(flair)
	flair = algorithms.loadMRI(flair)
	flair = algorithms.convertToArray(flair)
	print "\tFLAIR loading successful. " + str(time.clock() - start_time) + " seconds."
	
	#Get dimensions
	slices = len(flair)
	height = len(flair[0])
	width = len(flair[0][0])

	#Load T1 MRI (Conditional)
	t1 = []
	if len(dataset[3]) == 2 or len(dataset[3]) == 5:
		start_time = time.clock() 
		print "> Load T1 input."
		t1 = fileChooser('Load corresponding T1 input')
		t1 = algorithms.loadMRI(t1)
		t1 = algorithms.convertToArray(t1)
		print "\tT1 loading successful. " + str(time.clock() - start_time) + " seconds."

	#Load Brain Mask 
	start_time = time.clock()
	print "> Start brain mask creation."
	#mask = algorithms.createBrainMask(flair, range(startslice,endslice+1))
	mask = mask_const
	print "\tBrain mask creation successful. " + str(time.clock() - start_time) + " seconds."
	
	#KNN proper
	probabilities = knn(flair,t1,mask,dataset,algorithms)
	probabilities = probabilities[startslice:endslice+1]
	
	#Save KNN
	print "> Start saving KNN probabilistic slices."

	inputInfo = mricase + " [" + features + ", K_" + str(k) + "]"
	folder = "KNNResult\\" + inputInfo

	if os.path.exists(folder): shutil.rmtree(folder)
	try: os.makedirs(folder)
	except OSError as exc: raise

	image = sitk.GetImageFromArray(np.array(probabilities))
	filename = folder + "\\result.nhdr"
	sitk.WriteImage(image, filename)

	for i in range(len(probabilities)):
		imgplot = plt.imshow(probabilities[i], vmin=0, vmax=1)
		imgplot.axes.get_xaxis().set_visible(False)
		imgplot.axes.get_yaxis().set_visible(False)
		imgplot.set_cmap('spectral')
		plt.colorbar()

		filename = folder + "\\Slice" + str(i+startslice) + ".png"
		if os.path.isfile(filename): os.remove(filename)		
		plt.savefig(filename)
		plt.close()

	print "\tKNN probabilistic slices saved. ", time.clock() - start_time, "seconds."

	#Thresholding
	for i in range(len(mask)):
		for i1 in range(mask[i][1], mask[i][2]+1):
			for i2 in range(mask[i][3], mask[i][4]+1):
				if probabilities[i][i1][i2] < threshold1:
					probabilities[i][i1][i2] = 0.0

	#MRF
	global x, xa, xb, energy
	results = []
	start_time = time.clock()
	print "> Start MRF. [ from slice", startslice, "to", endslice, "]"

	#MAP iteration
	energy_arr = []
	for iteration in range(MAPiter):
		print "\tStart MAP iteration ", iteration+1

		#Slice iteration with brain mask
		energy = 0
		for i in range(len(probabilities)):
			x = probabilities[i]
			xa = None
			xb = None

			if i-1 >= 0: xa = probabilities[i-1]
			if i+1 < len(probabilities): xb = probabilities[i+1]

			print "\t\ Processing slice#", mask[i][0], "...\r",
			execute(mask[i][1], mask[i][2], mask[i][3], mask[i][4])
			probabilities[i] = x
	
		print "\t\ Processing slice#", endslice, "...done"
		print "\t\ Energy value: ", energy
		if len(energy_arr) != 0: print "\t\ Energy difference: ", energy_arr[-1] - energy
		else: print "\t\ Energy difference: n/a"
		energy_arr.append(energy)

	print "\tMRF successful. ", time.clock() - start_time, "seconds."
	
	global threshold2
	#Saving spectral
	start_time = time.clock()
	print "> Start saving probabilistic slices."

	inputInfo = mricase + " [" + features + ", K_" + str(k) + ", T1_" + str(threshold1) + ", I_" + str(iteration+1) + ", T2_" + str(threshold2) + "]"
	folder = "ProbabilisticMRI\\" + inputInfo

	if os.path.exists(folder): shutil.rmtree(folder)
	try: os.makedirs(folder)
	except OSError as exc: raise

	image = sitk.GetImageFromArray(np.array(probabilities))
	filename = folder + "\\result.nhdr"
	sitk.WriteImage(image, filename)

	for j in range(len(probabilities)):
		imgplot = plt.imshow(probabilities[j], vmin=0, vmax=1)
		imgplot.axes.get_xaxis().set_visible(False)
		imgplot.axes.get_yaxis().set_visible(False)
		imgplot.set_cmap('spectral')
		plt.colorbar()

		filename = folder + "\\Slice" + str(j+startslice) + ".png"
		if os.path.isfile(filename): os.remove(filename)		
		plt.savefig(filename)
		plt.close()

	print "\tProbabilistic slices saved. ", time.clock() - start_time, "seconds."

	while threshold2 < 1.01:
		#Thresholding
		for j in range(len(mask)):
			for i1 in range(mask[j][1], mask[j][2]+1):
				for i2 in range(mask[j][3], mask[j][4]+1):
					if probabilities[j][i1][i2] < threshold2:
						probabilities[j][i1][i2] = 0.0

		#Save segmented
		start_time = time.clock()
		print "> Start saving segmented slices."

		inputInfo = mricase + " [" + features + ", K_" + str(k) + ", T1_" + str(threshold1) + ", I_" + str(iteration+1) + ", T2_" + str(threshold2) + "]"
		folder = "SegmentedMRI\\" + inputInfo
		filename = folder + "\\result.nhdr"

		if os.path.exists(folder): shutil.rmtree(folder)
		try: os.makedirs(folder)
		except OSError as exc: raise

		image = sitk.GetImageFromArray(np.array(probabilities))
		sitk.WriteImage(image, filename)

		for j in range(len(probabilities)):
			filename = folder + "\\Slice" + str(j+startslice) + ".png"
			scipy.misc.imsave(filename, np.array(probabilities[j]))

		print "\tSegmented slices saved. ", time.clock() - start_time, "seconds."

		threshold2 = threshold2 + .01
	
	#Plot energy minimization
	print "> Plot energy minimization."
	plt.plot(range(1, len(energy_arr)+1), energy_arr, marker = 'o', linestyle = '--', color  ='r')

	folder = "EnergyMinimization\\"
	if not os.path.exists(folder):
		try: os.makedirs(folder)
		except OSError as exc: raise

	filename = folder + inputInfo + ".png"
	if os.path.isfile(filename): os.remove(filename)
	plt.savefig(filename)

	print ">> Done! Your MR image is segmented."


def main(params):
	global k, threshold1, MAPiter, threshold2

	if len(params) != 8: 
		print "Error number of arguments. Expected is 8."
		print "KNN_Test.py -k <# of neighbors> -t1 <threshold1> -i <MAP iterations> -t2 <threshold2>"
	
		return -1
		
	else: 
		for i in range(len(params)):
			try:
				if params[i] == '-k':
					k = int(params[i+1])
					if k <= 0:
						print "Invalid value of k. It should be greater than 0."
						return -1
				
				elif params[i] == '-t1' :
					threshold1 = float(params[i+1])
					if threshold1 > 1.0 or threshold1 < 0.0:
						print "Invalid value of threshold2. It should be between 0 and 1."
						return -1

				elif params[i] == '-i':
					MAPiter = int(params[i+1])
					if MAPiter <= 0:
						print "Invalid value of MAP iteration. It should be a positive integer."
						return -1

				elif params[i] == '-t2' :
					threshold2 = float(params[i+1])
					if threshold2 > 1.0 or threshold2 < 0.0:
						print "Invalid value of threshold2. It should be between 0 and 1."
						return -1

			except Exception as e: 	
				print "There's an error in parsing of the parameters. Review your input."
				print str(e)
				break 
	
	if k == -1 or threshold1 == -1 or MAPiter == -1 or threshold2 == -1: 
		print "Error. Missing parameter(s)."
		return -1 
	
	knn_mrf()

	
if __name__ == "__main__" : 
	main(sys.argv[1:])