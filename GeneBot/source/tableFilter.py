import argparse
import os.path

outputDirName = "output"

def makePath(fileName):
	return outputDirName + "/" + fileName


filterdTablePath = makePath("filtered_analysis_table.txt")


def filterTable(tableFile, filteredFile, scoreThreshold, igCount):
	
	filteredFile.write(tableFile.readline())

	predicates = [None, lambda count, score: count <= igCount,
						lambda count, score: score > scoreThreshold,
						lambda count, score: count <= igCount or score > scoreThreshold]
	predicate = predicates[(bool(scoreThreshold) << 1) + bool(igCount)]

	if scoreThreshold:
		scoreThreshold = str(scoreThreshold) + "%" # convert to percentage string

	lastSym = None
	guideCount = 0
	for (i, row) in enumerate(tableFile):
		row = row.rstrip()
		cells = row.split("\t")

		if cells[0] != lastSym:
			guideCount = 1
		else:
			guideCount += 1

		if predicate(guideCount, cells[1]):
			filteredFile.write(row + "\n") 
		
		lastSym = cells[0]
		print("\r%d guides have be processed" % i, end = "")

	print()

def doubleQuoteText(text):
	return '"' + text + '"'

def makeOutputFolder():
	if not os.path.exists(outputDirName):
		print('see the output files in the %s folder' % doubleQuoteText(outputDirName))
		os.mkdir(outputDirName)

def readOptions():
	des = """This is an program that filter the analyzed guides in the output table.
			If two options are specfied, it will enter the smart mode,
			which keeps the specfied number of guides first,
			and then continues fetching guides if their scores are above the threshold.
			Otherwise, it filters with single condition only."""
	parser = argparse.ArgumentParser(description = des, formatter_class = argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("filename", help = "specify a table file to process")
	parser.add_argument("-t", "--score-threshold", help = "specify the score threshold which filter the guides", type = int)
	parser.add_argument("-c", "--guide-count", help = "specify the ideal guides count that each symbol takes", type = int)

	args = parser.parse_args()

	return (args.filename, args.score_threshold, args.guide_count)


def main():
	(fileName, scoreThreshold, guideCount) = readOptions()

	print("Table Filter Starting...")
	makeOutputFolder()

	if not scoreThreshold and not guideCount:
		print("you should at least specify one option to process")
		return
			
	with open(fileName) as tableFile, open(filterdTablePath, "w") as filteredFile:
		filterTable(tableFile, filteredFile, scoreThreshold, guideCount);
		
	print("Process Complete. see output table at " + doubleQuoteText(filterdTablePath))
	
if __name__ == "__main__":
	main()