import SimpleITK as sitk
import sys 

class Algorithms: 

	def __init__(self):
	
		print "===Automatic Segmentation of Multiple Sclerosis Lesions using KNN and MRF==="

	#Loads the MRI to memory. 
	#Returns SITK Image of MRI. 
	def loadMRI(self,filename):

		sys.stdout.write("\tLoading MRI of " + filename + "...")
		
		image = sitk.ReadImage(filename)
		pixelID = image.GetPixelIDValue()
		caster = sitk.CastImageFilter()
		caster.SetOutputPixelType( pixelID )
		image = caster.Execute( image )
			
		print "done"
		
		return image
	
	#Returns numpy array of SITK image
	def convertToArray(self,image):
		
		return sitk.GetArrayFromImage(image)
	
	#Normalization method
	#Sets data to be in [0,1] scale
	def normalize(self,x,min,max):
	
		return (float(x)-float(min))/(float(max)-float(min))
	
	#Standardization method
	#Used to make mean be equal to 0 and a variance of 1
	def standardize(self,x,mean,std):
	
		return (float(x)-float(mean))/float(std)
		
	#Parser of the training data. 
	#will return an array A containing features,labels,technique,min/mean,max/std,
	# A[0] - features, A[1] - labels, A[2] - contains ['N'] if "normalized" or ['S'] if "standardized" , 
	# A[3] - min if A[2] is normalized else mean, A[4] - max if A[2] is normalized else std
	def readTrainingData(self,filename):
	
		sys.stdout.write("\tReading training data of " + filename + "...")


		file = open(filename,'r')
		
		#Features,5
		#Normalized,Standardized
		#Mean/Min
		#Max/Std
		
		featuresFound = False 
		techniqueFound = False 
		firstP_Found = False
		secondP_Found = False 
		normalize = False 
		standardize = False 
		
		featuresNum = 0
		technique = ''
		firstP_arr = []
		secondP_arr = []
		features = []
		labels = []
		
		error = False
		
		for line in file:
		
			data = line.split(",")
			
			if not featuresFound:
			
				if not( len(data) > 0 and len(data) < 3 ): 
				
					print "Error: Features line has invalid number of input."
					error = True 
					break 
			
				if data[0] == "#Features":
				
					featuresNum = int(data[1])
					
					if featuresNum <= 0:
					
						print "Error: Feature number is unacceptable."
						error = True 
						break 
					
					featuresFound = True
				else:
					
					print "Error: #Features line in 1st line is expected."
					error = True 
					break 
			
			elif not techniqueFound:
			
				if len(data) != 1: 
				
					print "Error: Invalid number of inputs in 2nd line."
					error = True 
					break 
			
				if data[0] == "#Normalized\n":
				
					normalize = True
					
				elif data[0] == "#Standardized\n":
				
					standardize = True 
					
				else :
				
					print "Error: Technique was expected on 2nd line. Input #Normalized or #Standardized"
					error = True 
					break 
			
				techniqueFound = True 
			
			elif not firstP_Found:
			
				if len(data[1:len(data)]) != featuresNum: 
				
					print "Error: Number of paramaters is not sufficient for 3rd line."
					error = True 
					break 
			
				if normalize: 
				
					if len(data[1:len(data)]) != featuresNum:
					
						print "Error: Number of parameters does not match number of features."
						error = True 
						break
				
					if data[0] == "#Min":
										
						for i in range(1,featuresNum+1): 
							
							firstP_arr.append(float(data[i]))
					
						firstP_Found = True 
					
					else: 
					
						print "Error: #Min was expected on 3rd line."
						error = True 
						break
				
				elif standardize:
				
					if len(data[1:len(data)]) != featuresNum:
					
						print "Error: Number of parameters does not match number of features."
						error = True 
						break
				
					if data[0] == "#Mean":
										
						for i in range(1,featuresNum+1): 
							
							firstP_arr.append(float(data[i]))
					
						firstP_Found = True 
					
					else: 
					
						print "Error: #Mean was expected on 3rd line."
						error = True 
						break 
				
			elif not secondP_Found:
			
				if len(data[1:len(data)]) != featuresNum: 
				
					print "Error: Number of paramaters is not sufficient for 4th line."
					error = True 
					break 
			
				if normalize: 
				
					if len(data[1:len(data)]) != featuresNum:
					
						print "Error: Number of parameters does not match number of features."
						error = True 
						break
				
					if data[0] == "#Max":
										
						for i in range(1,featuresNum+1): 
							
							secondP_arr.append(float(data[i]))
					
						secondP_Found = True 
					
					else: 
					
						print "Error: #Max was expected on 3rd line."
						error = True 
						break
				
				elif standardize:
				
					if len(data[1:len(data)]) != featuresNum:
					
						print "Error: Number of parameters does not match number of features."
						error = True 
						break
				
					if data[0] == "#Std":
										
						for i in range(1,featuresNum+1): 
							
							secondP_arr.append(float(data[i]))
					
						secondP_Found = True 
					
					else: 
					
						print "Error: #Std was expected on 3rd line."
						error = True 
						break 		

			elif featuresFound and techniqueFound and firstP_Found and secondP_Found: 

				if len(data) != featuresNum+1:
				
					print "Error: Number of features is lacking/excessive or label is not included."
					error = True 
					break 
				
				for i in range(len(data)-1):
				
					data[i] = float(data[i]) 
				
				data[len(data)-1] = int(data[len(data)-1])
				
				features.append(data[0:featuresNum])
				labels.append(data[featuresNum])
			
			else :
			
				error = True 
				break 
		
		if error :
	
			return -1 
		
		if normalize:
			technique = 'N'
		elif standardize:
			technique = 'S'
		
		print "done"
		print "\tTraining data successfully parsed."
		return [features,labels,technique,firstP_arr,secondP_arr]
	
	#The brain mask algorithm.
	#returns slice and corresponding rectangle coordinates of the brain mask.
	def createBrainMask(self,array, slice = [-1] ):

		sys.stdout.write("\tCreating brain mask...")
		mask = []
		
		if slice[0] == -1:
			length = range(len(array))
		else:
			length = slice
		
		for z in length:

			minY = -1
			maxY = -1
			minX = -1 
			maxX = -1
			
			#GET Min Y 
			for y in range(0,len(array[0])):
			
				for x in range(0,len(array[0][0])):
				
					if array[z][y][x] != 0:
					
						minY = y 
						break 
				
				if minY > -1 :
					
					break
			
			if minY != -1:
			
				#GET Max Y
				for y in range(minY,len(array[0])):
					
					background = True
					
					for x in range(0,len(array[0][0])):
					
						if array[z][y][x] != 0:
						
							background = False
							break 
					
					if background: 
						maxY = y-1 
						break 
						
				#GET Min X
				for x in range(0,len(array[0][0])):
				
					for y in range(minY,maxY+1):
					
						if array[z][y][x] != 0:
						
							minX = x
							break 
					
					if minX > -1 :
						break
				
				#GET Max X
				for x in range(minX,len(array[0][0])):
					
					background = True
					
					for y in range(minY,maxY+1):
					
						if array[z][y][x] != 0:
						
							background = False
							break 
					
					if background: 
						maxX = x-1 
						break 
				
				mask.append([z,minY,maxY,minX,maxX])
		
		print "done"			
		return mask 