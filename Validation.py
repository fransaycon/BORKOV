import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

class Validation:


	def __init___(self):
	
		print "Loaded Validation."
		
	def colorMap(self,arr):
	
		imgplot = plt.imshow(np.array(arr))
		imgplot.set_cmap('spectral')
		plt.colorbar()
		plt.show()
	
	def graphAccuracy(self,mask,arr,groundTruth,filename="Statistic Results",threshold,title):
	
		file = open(filename,'a')
		
		false_nega = 0
		false_pos = 0
		true_pos = 0
		true_nega = 0
		
		for i in range(len(mask)):
		
			for y in range(mask[i][1],mask[i][2]+1):
		
				for x in range(mask[i][3],mask[i][4]+1):
		
					if groundTruth[i][y][x] == 1 and arr[i][y][x] >= threshold:
						true_pos+=1 
					elif groundTruth[i][y][x] == 1 and arr[i][y][x] < threshold:
						false_nega+=1
					elif groundTruth[i][y][x] == 0 and arr[i][y][x] >= threshold:
						false_pos+=1 
					else: 
						true_nega+=1 
		
		total = false_nega+false_pos+true_nega+true_pos
		
		file.write(title+"\n")
		file.write("FALSE NEGATIVES: " + str(false_nega/total) + " TOTAL OF " + str(false_nega) + "\n")
		file.write("TRUE NEGATIVES: " + str(true_nega/total) + " TOTAL OF " + str(true_nega) + "\n")
		file.write("FALSE POSITIVES: " + str(false_pos/total) + " TOTAL OF " + str(false_pos) + "\n")
		file.write("TRUE POSITIVES " + str(true_pos/total) + " TOTAL OF " + str(true_pos) + "\n")