import SimpleITK as sitk
import numpy as np
import gc
import sys
import random
import math
from Algorithms import Algorithms

#this needs Settings.txt, loads the directories.
#Datasets_Raw_UNC: Index-0
#Datasets_Raw_CHB: Index-1
#Datasets_BrainMasks_UNC: Index-2
#Datasets_BrainMasks_CHB: Index-3
#Data_KNN_UNC: Index-4
#Data_KNN_CHB: Index-5 
def load():

	sys.stdout.write("Getting directories....")
	directories = [""] * 6

	try: 
		
		settings = open("Settings.txt",'r')
		
		for line in settings:
		
			data = line.split()
			
			if len(data) != 0 :
			
				if data[0][0] != '#' : 
				
					if data[0] == "Datasets_Raw_UNC:" :
						directories[0] = data[1]
					elif data[0] == "Datasets_Raw_CHB:" :
						directories[1] = data[1]
					elif data[0] == "Datasets_BrainMasks_UNC:" :
						directories[2] = data[1]
					elif data[0] == "Datasets_BrainMasks_CHB:" :
						directories[3] = data[1]
					elif data[0] == "Data_KNN_UNC:" :
						directories[4] = data[1]
					elif data[0] == "Data_KNN_CHB:" :
						directories[5] = data[1]
						
	except Exception, e: 
		print e
	
	print "done"
	
	return directories

def createMask(data,filename):

	sys.stdout.write("Creating brain mask of " + filename + ".....")

	text_file = open( filename, "w")
	voxelNum = 0
	
	for x in range(len(data)):
		for y in range(len(data[0])):
			for z in range(len(data[0][0])):
			
				if data[x][y][z] != 0 :
					text_file.write(str(x) + " " + str(y) + " " + str(z) + "\n" )
					voxelNum+=1
	
	print "done"
	print "Total Tissue Voxles: " + str(voxelNum)
	
	text_file.close()

def separate(groundTruth,filename,f1,f2):

	sys.stdout.write("Separation process.....")

	sep = open(filename,'r')
	lesion = open(f1,'w')
	nonLesion = open(f2,'w')
	
	numL = 0.0
	numNL = 0.0
	
	for line in sep:
	
		data = line.split()
		x = int(data[0]) 
		y = int(data[1])
		z = int(data[2])
		
		if(groundTruth[x][y][z] == 1):
			lesion.write( str(x) + " " + str(y) + " " + str(z) + "\n" )
			numL += 1.0
		else:
			nonLesion.write( str(x) + " " + str(y) + " " + str(z) + "\n" )
			numNL += 1.0
			
	sep.close()
	lesion.close()
	nonLesion.close()
	
	print "done"
	
	return (numL/numNL)

def getSamples(f1,f2,percent):

	sys.stdout.write("Getting samples per nonlesion set of " + f1 + ".....")

	population = open(f1,'r')
	sample = open(f2,'w')
	
	temp = []
	iter = 0
	
	for line in population:
	
		data = line.split()
		temp.append(data)
		iter+=1
		
		if iter == 500000:
		
			temp = random.sample(temp,int(percent*iter))
			
			for i in range(len(temp)):
			
				sample.write(temp[i][0] + " " + temp[i][1] + " " + temp[i][2] + "\n")
			
			temp = []
			iter = 0
			gc.collect()
	
	if iter > 0 :
	
		temp = random.sample(temp,int(percent*iter))
				
		for i in range(len(temp)):
		
			sample.write(temp[i][0] + " " + temp[i][1] + " " + temp[i][2] + "\n")
		
		temp = []
		iter = 0
	
	population.close()
	sample.close()
		
	print "done"

def createChunks(flair,t1,lesion_f,nonLesion_f,filename,percentage,algorithms,option):

	sys.stdout.write("Creating chunk of " + filename + "....." )

	lesion = open(lesion_f,'r')
	nonLesion = open(nonLesion_f,'r')
	file = open(filename,'w')
	
	lesion_sample = []
	nonLesion_sample = []
	
	for line in lesion:
		
		data = line.split()
		lesion_sample.append(data)
	
	lesion_sample = random.sample(lesion_sample,int(percentage*len(lesion_sample)))
	
	for data in lesion_sample:
	
		x = int(data[0])
		y = int(data[1]) 
		z = int(data[2]) 
		
		if option==1:
			file.write(str(flair[x][y][z]) + ",1\n")
		elif option == 2:
			file.write(str(flair[x][y][z]) + "," + str(t1[x][y][z])+ ",1\n")	
		elif option == 3: 
			file.write(str(x) + "," + str(y) + "," + str(z) + "," + str(flair[x][y][z]) + ",1\n")
		elif option == 4:
			file.write(str(x) + "," + str(y) + "," + str(z) + "," + str(flair[x][y][z]) + "," + str(t1[x][y][z])+ ",1\n")
		
	lesion.close()
	lesion_sample = []
	gc.collect()
	
	for line in nonLesion:
		
		data = line.split()
		nonLesion_sample.append(data)
	
	'''numpy can only handle much data. in the getSample method we get random samples every 500000. to the percentage of lesion and nonlesion per dataset
	this might induce a difference in the population of nonlesion and lesion when it comes to the samples.
	wherein the nonlesions are behind. so we increase a the percentage to get in every chunck by some amount to compensate .'''
	nonLesion_sample = random.sample(nonLesion_sample,int((percentage+0.008)*len(nonLesion_sample)))
	
	for data in nonLesion_sample:

		x = int(data[0])
		y = int(data[1]) 
		z = int(data[2]) 
		
		if option==1:
			file.write(str(flair[x][y][z]) + ",0\n")
		elif option == 2:
			file.write(str(flair[x][y][z]) + "," + str(t1[x][y][z])+ ",0\n")	
		elif option == 3: 
			file.write(str(x) + "," + str(y) + "," + str(z) + "," + str(flair[x][y][z]) + ",0\n")
		elif option == 4:
			file.write(str(x) + "," + str(y) + "," + str(z) + "," + str(flair[x][y][z]) + "," + str(t1[x][y][z])+ ",0\n")
	
	nonLesion.close()
	nonLesion_sample = []
	gc.collect()
	
	print "done"
	

#combine all chunks 
def createKnowledge(filename,dir):

	sys.stdout.write("Combining all chunks....")

	brain = open(filename,'w')
	
	for i in range(7):
	
		chunk = open(dir[4] + "\\UNC_CASE0" + str(i+1) + "_CHUNK.txt",'r')
		
		for line in chunk:
		
			brain.write(line)
	
	for i in range(7):
	
		chunk = open(dir[5] + "\\CHB_CASE0" + str(i+1) + "_CHUNK.txt",'r')
		
		for line in chunk:
		
			brain.write(line)
	
	brain.close()
	
	print("done")

def main(params):

	option = 0 
	samplePercentage = 0.0
	
	if len(params) != 4 :
		print "Error: expected length of parameters is 4."
		print "KNN_Train.py -o [1,4] -p percent"
		print "1 - FLAIR Intensity Only"
		print "2 - FLAIR and T1 Intensity"
		print "3 - x,y,z and FLAIR Intensity "
		print "4 - x,y,z, FLAIR and T1 Intensity"
	else: 
		for i in range(len(params)):
		
			try:
		
				if params[i] == '-o':
				
					option = int(params[i+1])
					
					if option > 4 or option < 1:
					
						print "Invalid value of option. Choose between 1 and 4"
						return -1
				
				elif params[i] == '-p' :
				
					samplePercentage = float(params[i+1])
					
					if samplePercentage > 1.0 or samplePercentage < 0:
					
						print "Invalid value of percentage. Percentage should be between 0 and 1."
						return -1
			
			except Exception as e: 
				
				print "There's an error in parsing of the parameters. Review your input."
				print str(e)
				break 
	
	print ("Program has started.")
	#load all file dependencies
	
	directories = load()
	
	sys.stdout.write("Producing filenames.....")
	brainmasks_UNC = ["UNC_CASE0" + str(n) + "_BRAINMASK.nhdr" for n in range(1,8)]
	brainmasks_CHB = ["CHB_CASE0" + str(n) + "_BRAINMASK.nhdr" for n in range(1,8)]
	groundTruth_UNC = ["UNC_train_Case0" + str(n) + "_lesion.nhdr" for n in range(1,8)]
	groundTruth_CHB = ["CHB_train_Case0" + str(n) + "_lesion.nhdr" for n in range(1,8)]
	flair_UNC = ["UNC_train_Case0" + str(n) + "_FLAIR.nhdr" for n in range(1,8)]
	flair_CHB = ["CHB_train_Case0" + str(n) + "_FLAIR.nhdr" for n in range(1,8)]
	t1_UNC = ["UNC_train_Case0" + str(n) + "_T1.nhdr" for n in range(1,8)]
	t1_CHB = ["CHB_train_Case0" + str(n) + "_T1.nhdr" for n in range(1,8)]
	print "done"
	sys.stdout.write("Initializing subdirectories.....")
	subdirs_UNC = ["UNC_CASE0" + str(n) for n in range(1,8)]
	subdirs_CHB = ["CHB_CASE0" + str(n) for n in range(1,8)]
	print "done"
	
	percentages = []
	algorithms = Algorithms()
	
	'''
	#create the 5 masks in UNC 
	for i in range(7):
	
		data = algorithms.loadMRI(directories[2] + "\\" + subdirs_UNC[i] + "\\" + brainmasks_UNC[i])
		data = sitk.GetArrayFromImage(data)
		createMask(data,directories[4] + "\\UNC_CASE0" + str(i+1) + "_BRAINMASK.txt")
		data = []
		gc.collect() #flush memory
	
	#create the 5 masks in CHB 
	for i in range(7):
	
		data = algorithms.loadMRI(directories[3] + "\\" + subdirs_CHB[i] + "\\" + brainmasks_CHB[i])
		data = sitk.GetArrayFromImage(data)
		createMask(data,directories[5] + "\\CHB_CASE0" + str(i+1) + "_BRAINMASK.txt")
		data = []
		gc.collect() #flush memory
	
	
	#separate into lesions and nonlesions UNC 
	for i in range(7):
	
		groundTruth = algorithms.loadMRI(directories[0] + "\\" + subdirs_UNC[i] + "\\" + groundTruth_UNC[i] )
		groundTruth = sitk.GetArrayFromImage(groundTruth)
		percentages.append(separate(groundTruth,directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_BRAINMASK.txt",directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_LESIONS.txt",directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_NONLESIONS.txt"))
		groundTruth = []
		gc.collect() #flush memory 
	
	#separate into lesions and nonlesions CHB
	for i in range(7):
	
		groundTruth = algorithms.loadMRI(directories[1] + "\\" + subdirs_CHB[i] + "\\" + groundTruth_CHB[i] )
		groundTruth = sitk.GetArrayFromImage(groundTruth)
		percentages.append(separate(groundTruth,directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_BRAINMASK.txt",directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_LESIONS.txt",directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_NONLESIONS.txt"))
		groundTruth = []
		gc.collect() #flush memory 
	
	#get samples from nonlesion UNC
	for i in range(7):
	
		getSamples(directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_NONLESIONS.txt",directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_SAMPLES.txt", percentages[i])
	
	#get samples from nonlesion CHB
	for i in range(7):
	
		getSamples(directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_NONLESIONS.txt",directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_SAMPLES.txt", percentages[i+5])

	'''
	
	#create the chunks of knowledge per lesion and sample UNC
	for i in range(7):
		
		flair = algorithms.loadMRI(directories[0] + "\\" + subdirs_UNC[i] + "\\" + flair_UNC[i] )
		flair = sitk.GetArrayFromImage(flair)
		t1 = []
		
		if option == 2 or option == 4 :
			t1 = algorithms.loadMRI(directories[0] + "\\" + subdirs_UNC[i] + "\\" + t1_UNC[i] )
			t1 = sitk.GetArrayFromImage(t1)
			
		createChunks(flair,t1,directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_LESIONS.txt",directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_SAMPLES.txt",directories[4] + "\\" + "UNC_CASE0" + str(i+1) + "_CHUNK.txt",samplePercentage,algorithms,option)
		flair = []
		t1 = []
	
	#create the chunks of knowledge per lesion and sample CHB
	for i in range(7):
		
		flair = algorithms.loadMRI(directories[1] + "\\" + subdirs_CHB[i] + "\\" + flair_CHB[i] )
		flair = sitk.GetArrayFromImage(flair)
		t1 = [] 
		
		if option == 2 or option == 4 :
		
			t1 = algorithms.loadMRI(directories[1] + "\\" + subdirs_CHB[i] + "\\" + t1_CHB[i] )
			t1 = sitk.GetArrayFromImage(t1)
			
		createChunks(flair,t1,directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_LESIONS.txt",directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_SAMPLES.txt",directories[5] + "\\" + "CHB_CASE0" + str(i+1) + "_CHUNK.txt",samplePercentage,algorithms,option)
		flair = []
		t1 = []
	
	#combine all chunks 
	createKnowledge("Brain.csv",directories)
	
if __name__ == "__main__" : 
	main(sys.argv[1:])