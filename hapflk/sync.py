import collections
import numpy as np
import gzip

### CLASS DEFINITIONS ###

class SinglePop:
	def __init__(self,a,t,c,g,n,deletion):
		"""
		Represents the allele counts of a single population
		"""
		self.__a=int(a)
		self.__t=int(t)
		self.__c=int(c)
		self.__g=int(g)
		self.__n=int(n)
		self.__del=int(deletion)
		
	# changed to return true or false
	def issnp(self,mincount):
		"""
		>>> SinglePop(2,0,0,2,0,0).issnp(1)
		True
		>>> SinglePop(2,0,0,2,0,0).issnp(3)
		False
		>>> SinglePop(2,0,0,1,0,0).issnp(2)
		False
		>>> SinglePop(2,0,0,0,0,0).issnp(1)
		False
		>>> SinglePop(2,2,0,0,0,0).issnp(1)
		False
		"""
		alcount=self.count_alleles(mincount)
		if(alcount>1):
			return True
		else:
			return False

	def count_alleles(self,mincount):
		"""
		>>> SinglePop(2,0,0,2,0,0).count_alleles(1)
		2
		>>> SinglePop(2,0,0,2,0,0).count_alleles(2)
		2
		>>> SinglePop(2,0,0,2,0,0).count_alleles(3)
		0
		"""
		alcount=0
		if self.A>=mincount:
			alcount+=1
		if self.T>=mincount:
			alcount+=1
		if self.C>=mincount:
			alcount+=1
		if self.G>=mincount:
			alcount+=1
		return alcount

	def countForAllele(self,allele):
		"""
		>>> SinglePop(2,0,0,0,0,0).countForAllele("A")
		2
		>>> SinglePop(2,4,0,0,0,0).countForAllele("T")
		4
		"""
		return eval("self."+allele)		

	@property
	def A(self):
		"""
		>>> SinglePop(2,0,0,0,0,0).A
		2
		"""
		return self.__a
	
	@property
	def Af(self):
		"""
		>>> SinglePop(4,0,0,0,0,0).Af
		1.0
		>>> SinglePop(4,4,0,0,0,0).Af
		0.5
		"""
		return float(self.A)/float(self.cov)
	
	
	@property
	def T(self):
		"""
		>>> SinglePop(2,0,0,0,0,0).T
		0
		>>> SinglePop(0,2,0,0,0,0).T
		2
		"""
		return self.__t
	
	@property
	def Tf(self):
		"""
		>>> SinglePop(0,4,0,0,0,0).Tf
		1.0
		>>> SinglePop(4,4,0,0,0,0).Tf
		0.5
		"""
		return float(self.T)/float(self.cov)
		

	@property
	def C(self):
		"""
		>>> SinglePop(0,0,2,0,0,0).C
		2
		"""
		return self.__c
	
	@property
	def Cf(self):
		"""
		>>> SinglePop(0,0,4,0,0,0).Cf
		1.0
		>>> SinglePop(0,4,4,0,0,0).Cf
		0.5
		"""
		return float(self.C)/float(self.cov)

	@property
	def G(self):
		"""
		>>> SinglePop(0,0,0,3,0,0).G
		3
		"""
		return self.__g
	
	@property
	def Gf(self):
		"""
		>>> SinglePop(0,0,0,4,0,0).Gf
		1.0
		>>> SinglePop(0,4,0,4,0,0).Gf
		0.5
		"""
		return float(self.G)/float(self.cov)

	@property
	def N(self):
		"""
		>>> SinglePop(0,0,0,0,5,0).N
		5
		"""
		return self.__n
		
	
	@property
	def deletion(self):
		"""
		>>> SinglePop(0,0,0,0,0,6).deletion
		6
		"""
		return self.__del

	@property
	def cov(self):
		"""
		>>> SinglePop(2,3,1,1,0,0).cov
		7
		"""
		return self.A+self.T+self.C+self.G

	@property
	def totcov(self):
		"""
		>>> SinglePop(2,3,1,1,2,1).totcov
		10
		"""
		return self.cov+self.N+self.deletion
	
	def __str__(self):
		"""
		>>> str(SinglePop(6,5,4,3,2,1))
		'6:5:4:3:2:1'
		"""
		return ":".join(map(str,[self.A,self.T,self.C,self.G,self.N,self.deletion]))



class PopLine:
	def __init__(self,chr,pos,refc,populations,pvalue=None):
		self.__chr=chr
		self.__pos=int(pos)
		self.__refc=refc
		self.__populations=populations
		self.__pvalue=pvalue
		

	@property
	def chr(self):
		"""
		>>> PopLine("2L",1,"N",[],0.2).chr
		'2L'
		"""
		return self.__chr
	
	@property
	def pos(self):
		return self.__pos
	
	@property
	def refc(self):
		return self.__refc

	@property
	def populations(self):
		return self.__populations
	
	def subpopulations(self,populations):
		tpops=self.populations
		toret=[]
		for i in populations:
			toret.append(tpops[i-1])
		return toret
	
	@property
	def popcount(self):
		return len(self.__populations)
	
	def __str__(self):
		"""
		"""
		popstr="\t".join(map(str,self.populations))
		tojoin=[self.chr,self.pos,self.refc,popstr]
		return "\t".join([str(x) for x in tojoin])
		

class SyncReader:
	def __init__(self,fileName):
		self.__filename=fileName
		if(isinstance(fileName,str)):
			if isGZIP(fileName):
				self.__filehandle=gzip.open(fileName,"r")
			else:
				self.__filehandle=open(fileName,"r")
		else:
			self.__filehandle=fileName

	def __iter__(self):
		return self
	
	def next(self):
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line != "":
				break
		
		a=line.split()
		chr=a.pop(0)
		pos=a.pop(0)
		refc=a.pop(0)
		population=[]
		for p in a:
			po=None
			if p=="-":
				po=SinglePop(0,0,0,0,0,0)
			else:
				s=p.split(":")
				po=SinglePop(*s)
			population.append(po)
		
		return PopLine(chr,pos,refc,population)

	def close(self):
		self.__filehandle.close()

### METHODS ###

def isGZIP(fileName):
	"""
	Check if a file ends with a gzip/bgzip extension
	"""
	return fileName.endswith("gz") or fileName.endswith("bgz")

def popCount(fileName):
	"""
	Compute the number of populations in a sync file
	"""
	syncReader = SyncReader(fileName)
	npops = syncReader.next().popcount
	syncReader.close()
	return npops

def commaListToPopulationIndexes(popList, maxPop):
	"""
	get the 1-based list of populations from a comma separated format

	popList: comma-separated string with indexes
	"""
	toret = list(sorted(set(map(int, popList.split(",")))))
	if 0 in toret:
		print "Population list should be 1-based and found a 0"
		raise
	elif toret[-1] > maxPop:
		print "Population list should could not contain more than the number of populations in the file"
		raise
	return toret
