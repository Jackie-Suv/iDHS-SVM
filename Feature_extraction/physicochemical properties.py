def detectingPythonVersionAndOutputHelpInfo():
	if sys.version[0] == '2':
		print("""\nVersionError: The current python version: '%s',
		You must use python version 3 or later to run this script!\n"""%((sys.version).split('\n')[0]))
		exit(0)
	else:
		pass

	try:
		if sys.argv[1] == "--help":
			printHelpInfo()
		else:
			pass
	except:
		printHelpInfo()

def obtainExternalParameters():
	try:
		if sys.argv[1] == "-t":
			typeOfPseKNC = int(sys.argv[2])
		else:
			assert 0
		if sys.argv[3] == "-w":
			weightFactor = float(sys.argv[4])
		else:
			assert 0
		if sys.argv[5] == "-r":
			lambdaPara = int(sys.argv[6])
		else:
			assert 0
		if sys.argv[7] == "-i":
			in_filename = sys.argv[8]
		else:
			assert 0
		if sys.argv[9] == "-o":
			out_filename = sys.argv[10]
		else:
			assert 0
	except:
			printHelpInfo()
	return typeOfPseKNC, weightFactor, lambdaPara, in_filename, out_filename


def obtainNucleotidesPhysicoChemicalDict(filename):
	phychemdict = dict()
	phyChemList = []

	count_line = 0
	f = open(filename)
	for eachline in f:
		count_line += 1
		temp = eachline.strip().split("\t")
		if count_line == 1:
			nucleotides = temp[1:]
		else:
			phyChemList.append(temp[0])
			phychemdict[temp[0]] = dict()
			count_num = 0
			for each in temp[1:]:
				count_num += 1
				phychemdict[temp[0]][nucleotides[count_num-1]] = float(each)
	f.close()

	return phychemdict, nucleotides, phyChemList


def generateCsvFormatLinebyType1PseKNC(in_file, out_file, numConjoin, weightFact, lambdaPara, nucleoStandDict, nucleotides, phyChemList):
	pass


def calculateFeatureValueByCorrFactorsDictAndOccurfrequencyType2(corrFactorsDict, occurfrequency, weightFact, lambdaPara, nucleotides, phyChemList):
	featureValueStr_kmer = ''
	featureValueStr_phy = ''

	corrSum = 0
	for eachkey in phyChemList:
		corrSum += sum(corrFactorsDict[eachkey].values())
	corrPart = corrSum*weightFact
	
	for eachNuc in nucleotides:
		featureValueStr_kmer += ',%.8f'%(occurfrequency[eachNuc])

	for i in range(lambdaPara):
		for phyChemKey in phyChemList:
			featureValueStr_phy += ',%.8f'%(weightFact*corrFactorsDict[phyChemKey][i+1])

	return featureValueStr_kmer,featureValueStr_phy

def calculateOccurenceFrequencyOfOlinucletide(kTuples, nucleotides):
	occurfrequency = dict()
	tupleLen = len(kTuples)
	for each in nucleotides:
		occurfrequency[each] = kTuples.count(each)/tupleLen
	
	return occurfrequency

def calculateAllCorrelationFactorAndOccurenceFrequencyType2(sequence, numConjoin, nucleoStandDict, nucleotides, phyChemList):
	corrFactorsDict = dict()

	seqLen = len(sequence)
	kTuples = []  # All twins in a sequence
	for i in range(seqLen-numConjoin+1):
		kTuples.append(sequence[i:i+numConjoin])

	occurfrequency = calculateOccurenceFrequencyOfOlinucletide(kTuples, nucleotides)


	for eachName in phyChemList:
		corrFactorsDict[eachName] = dict()
		for lamPa in range(1,(lambdaPara+1)):
			temp = []
			for kTuplesIndex in range(len(kTuples)-lamPa):
				preKTuple = kTuples[kTuplesIndex]
				backKTuple = kTuples[kTuplesIndex+lamPa]
				try:
					tempNumber = nucleoStandDict[eachName][preKTuple]*nucleoStandDict[eachName][backKTuple]
					temp.append(tempNumber)
				except:
					continue

			corrFactorsDict[eachName][lamPa] = sum(temp)/len(temp)

	return corrFactorsDict, occurfrequency


def generateCsvFormatNoteLineType2_kmer(nucleotides):
	noteLine = 'class'
	for eachNuc in nucleotides:
		noteLine += ",%s_f"%(eachNuc)

	return noteLine+'\n'


def generateCsvFormatNoteLineType2_phy(lambdaPara, phyChemList):
	noteLine = 'class'

	for i in range(lambdaPara):
		for eachName in phyChemList:
			noteLine += ",%s_%d%d"%(eachName,(i+1),lambdaPara)

	return noteLine+'\n'


def generateCsvFormatLinebyType2PseKNC(in_file, out_file, numConjoin, weightFact, lambdaPara, nucleoStandDict, nucleotides, phyChemList):
	basename, address = os.path.splitext(out_file)
	out_file_kmer = basename + '_kmer' + address
	out_file_phy = basename + '_phy' + address

	g1 = open(out_file_kmer,'w')
	g2 = open(out_file_phy,'w')
	g1.write(generateCsvFormatNoteLineType2_kmer(nucleotides))
	g2.write(generateCsvFormatNoteLineType2_phy(lambdaPara, phyChemList))
	g1.close()
	g2.close()


	g1 = open(out_file_kmer,'a')
	g2 = open(out_file_phy,'a')
	f = open(in_file)
	count_line = 0
	for eachline in f:
		if eachline[0] == '>':
			if count_line in list(range(0,280)):
				sampleType = 1
			else:
				sampleType = 2
			count_line += 1
			#sampleType = re.findall(r'>(\d)', eachline)[0] #re.findall(r'@(\d)@', eachline)[0]
		else:
			sequence = eachline.strip()
			[corrFactorsDict, occurfrequency] = calculateAllCorrelationFactorAndOccurenceFrequencyType2(sequence, numConjoin, nucleoStandDict, nucleotides, phyChemList)
			[tempLine_kmer,tempLine_phy] = calculateFeatureValueByCorrFactorsDictAndOccurfrequencyType2(corrFactorsDict, occurfrequency, weightFact, lambdaPara, nucleotides, phyChemList)
			g1.write("%d%s\n"%(sampleType,tempLine_kmer))
			g2.write("%d%s\n"%(sampleType,tempLine_phy))
	g1.close()
	g2.close()
	f.close()

import re
import sys, os

detectingPythonVersionAndOutputHelpInfo()
[typeOfPse, weightFact, lambdaPara, in_file, out_file] = obtainExternalParameters()
numConjoin = 2 #2联体
if __name__ == '__main__':
	diNcleoStandFile = r'6_standard.txt'
	[nucleoStandDict, nucleotides, phyChemList] = obtainNucleotidesPhysicoChemicalDict(diNcleoStandFile)
	
	if typeOfPse == 1:
		generateCsvFormatLinebyType1PseKNC(in_file, out_file, numConjoin, weightFact, lambdaPara, nucleoStandDict, nucleotides, phyChemList)
	elif typeOfPse == 2:
		generateCsvFormatLinebyType2PseKNC(in_file, out_file, numConjoin, weightFact, lambdaPara, nucleoStandDict, nucleotides, phyChemList)
	else:
		printHelpInfo()

	print("------Finished!------")
