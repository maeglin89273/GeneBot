from urllib.parse import urlencode
from urllib.request import urlopen
from queue import Queue
from queue import Empty
from sys import stdout

import urllib
import argparse
import re
import sys
import os
import os.path
import json
import threading
import time
import operator

searchHost = "http://www.ncbi.nlm.nih.gov"
startPoint = "http://www.ncbi.nlm.nih.gov/gene/"
fastaUrl = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&val={0}&from={1}&to={2}"

crisprHost = "http://crispr.mit.edu"
crisprPostUrl = "/j/post_new"
analysisProfileUrl = crisprHost + "/job/"
gbDownloadUrl = "/export/guides_gb/{0}"

fastaLinkRegex = re.compile(r'<a title="Nucleotide FASTA report" href="(.+)" ref=".+">FASTA</a>')
multiResultHandleLinkRegex = re.compile(r'<a href="(/gene/\d+)" ref="ordinalpos=1')
uidRegex = re.compile(r'<meta name="ncbi_uidlist" content="(\d+)" />')
fromRegex = re.compile(r'/nuccore/.+?report=fasta&from=(\d+)&to=\d+')
symbolRegex = re.compile(r"([\w\-]+)\s+.+")
tableDataRegex = re.compile(r'protein_bind\s+(.+?)\n\s+/bound_moiety="CRISPR\s([\w\-]+?)\sguide\s\#\d+\s+([A-Z]+)".+?"score":\s"(\d+)%"', re.DOTALL)


outputDirName = "output"

def makePath(fileName):
	return outputDirName + "/" + fileName

logFilePath = makePath("Gene_Bot.log")
tableFilePath = makePath("analysis_table.txt")

total = 0
completeCount = 0
errorCount = 0
warningCount = 0

taskQueue = Queue()
responseQueue = Queue()

progressTextTmp = "\rdownloaded %d analyses"

log = None
INFO = "INFO"
WARNING = "WARNING"
ERROR = "ERROR"
logFileLock = threading.Lock()
logTemplate = "{} {}:\n\tSearch word:\n\t\t{}\n\t{} description:\n\t\t{}\n\n"

analysisTableFile = None
tableFileLock = threading.Lock()
rowTemplate =  "{}" + "\t{}" * 5 + "\n"

backupFile = None
backupFileLock = threading.Lock()

rangeOffset = 0
seqLength = 0
email = ""
watchingInterval = 0
watchingTimeout = 0
logVerbose = False

def storeTable(table):
	with tableFileLock:
		for row in table:
			analysisTableFile.write(rowTemplate.format( row["symbol"],
														str(row["score"]) + "%",
														row["pool"],
														row["sequence"],
														row["direction"],
														row["protein_bind"]))

		analysisTableFile.flush()

def parseSymbol(searchWord):
	return symbolRegex.search(searchWord).group(1)

def parseAnalysis(searchWord, analysis):
	table = []
	symbol = parseSymbol(searchWord)
	for match in tableDataRegex.finditer(analysis):
		row = { "symbol": symbol,
		 		"protein_bind": match.group(1),
		 		"direction": match.group(2),
		 		"sequence": match.group(3),
		 		"score": int(match.group(4)),
		 		"pool": "NUCLEAR"}

		table.append(row)

	return table

def storeAnalysis(searchWord, response):
	
	table = parseAnalysis(searchWord, getContentAndClose(response))
	table.sort(key = operator.itemgetter("score"), reverse = True)
	storeTable(table)

	finishTask()

	logProgress(INFO , searchWord, "analysis stored")

def watchAnalysis(searchWord, jobKey):
	resultUrl = crisprHost + gbDownloadUrl.format(jobKey)
	statusCode = 500
	errorMsg = ""
	response = None
	timeoutCount = watchingTimeout * 60 / watchingInterval if watchingTimeout else float("inf")

	logProgress(INFO, searchWord, "watching analysis progress")
	
	i = 0
	while statusCode == 500 and i < timeoutCount:
		time.sleep(watchingInterval)
		i += 1
		try:
			response = urlopen(resultUrl)
			statusCode = 200 #sever process complete
			logProgress(INFO, searchWord, "finish analyzing")
			break

		except Exception as e:
			statusCode = e.code if type(e) == urllib.error.HTTPError else -1
			errorMsg = str(e.reason) if type(e) == urllib.error.URLError else str(e)
	

	else: # non 500 error or waiting timeout, it have to be handled internally because it isn't invoked by a worker
		backupAnalysisProcess(searchWord, jobKey)
		if statusCode != 500:
			handleError(Exception("CRISPR server analyzing error occurs:\n" +
							 		errorMsg + ", code: " + str(statusCode) +
							 		"\nvisit " + analysisProfileUrl + jobKey + " for more details,\n" +
							 		"or run the backup shell which helps you to fetch the lost analysis"), searchWord)
		else: # waiting timeout
			errorMsg = "this worker spends too much time for waiting this analysis"
			handleError(Exception("Waiting timeout:\n" +
							 		errorMsg + ", code: " + str(-500) +
							 		"\nvisit " + analysisProfileUrl + jobKey + " for more details,\n" +
							 		"or run the backup shell which helps you to fetch the lost analysis"), searchWord)


	#submit new analysis request first, analysis response will be process later
	responseQueue.put((searchWord, response))


def backupAnalysisProcess(searchWord, jobKey):
	with backupFileLock:
		global backupFile
		if not backupFile:
			backupFile = open(makePath("analysis_backup.txt"), "w")

		backupFile.write(parseSymbol(searchWord) + " " + jobKey + "\n")
		backupFile.flush()

def closeBackupFile():
	if backupFile:
		backupFile.close()	


def submitToCRISPR(searchWord, sequence):
	payload = { "name": searchWord,
				"email": email,
				"inputRadios": "other_region",
				"genome": "hg19",
				"query": sequence}
	
	response = urlopen(crisprHost + crisprPostUrl, data = urlencode(payload).encode())
	status = json.loads(response.read().decode())
	if status["status"] != "success":
		raise Exception("CRISPR Server Error")

	logProgress(INFO, searchWord, "submitted to CRISPR server, visit:\n" + analysisProfileUrl + status["job_key"])

	#start watching the analysis progress
	t = threading.Thread(target = watchAnalysis, args = (searchWord, status["job_key"]))
	t.start()


def fetchSequence(searchWord, uidValue, fromValue, toValue):

	response = urlopen(fastaUrl.format(uidValue, fromValue, toValue))
	response.readline() #skip the header line
	sequence = response.read().decode()
	response.close()

	logProgress(INFO, searchWord, "sequence fetched")
	return (searchWord, sequence.rstrip())
	
	

def enterFastaLink(searchWord, link):

	link = link.replace('amp;', '')
	response = urlopen(searchHost + link)

	html = getContentAndClose(response, 1500)

	m = fromRegex.search(link)

	uidValue = uidRegex.search(html).group(1)
	fromValue = int(m.group(1)) + rangeOffset
	toValue = fromValue + seqLength - 1

	return fetchSequence(searchWord, uidValue, fromValue, toValue)

def enterFirstSearchResult(searchWord, html):
	global warningCount
	warningCount = warningCount + 1

	logProgress(WARNING, searchWord, "Multiple search results, we choose the first one")

	firstResultLink = multiResultHandleLinkRegex.search(html).group(1)
	response = urlopen(searchHost + firstResultLink)
	return getContentAndClose(response) 

def searchKeyword(searchWord):

	logProgress(INFO, searchWord, "start searching")

	response = urlopen(startPoint + "?" + urlencode({"term": searchWord}))
	html = getContentAndClose(response)

	m = fastaLinkRegex.search(html)

	if not m:
		m = fastaLinkRegex.search(enterFirstSearchResult(searchWord, html))
		
	return enterFastaLink(searchWord, m.group(1))

		

def getContentAndClose(response, amt = None):
	content = response.read(amt).decode()
	response.close()
	return content

def finishTask():
	global completeCount
	completeCount = completeCount + 1
	updateStatus()

def handleError(e, searchWord):
	global errorCount
	errorCount = errorCount + 1

	logProgress(ERROR, searchWord, str(e))

	updateStatus()


def updateStatus():
	progressText = progressTextTmp % completeCount
	errorText = ""
	percentageText = ""
	
	if errorCount > 0:
		errorText = ", skipped %d error(s)" % errorCount

	if total > 0:
		percentageText = " (%.2f%%)" % ((completeCount + errorCount) * 100 / total)

	progressText += percentageText + errorText

	stdout.write(progressText)
	stdout.flush()

def logProgress(type, searchWord, content):
	if not logVerbose and type == INFO:
		return

	with logFileLock:
		content = content.replace("\n", "\n\t\t")
		log.write(logTemplate.format(getTime(),type, searchWord, type, content))
		log.flush()

def hasTasks():
	return not (total and taskQueue.empty())
# this function shouldn't be modified because it optimizes the workflow that takes the shortest time
def worker():

	task = None
	lastSubmitted = False
	newSubmit = False # it will be set to lastSubmitted later
	while hasTasks() or lastSubmitted:
		submitArgs = None
		# stop searching if every keyword is searched
		# it prevent the lock caused by taskQueue.get()
		if hasTasks():
			try:
				task = taskQueue.get(timeout = 30) # start the next workflow. If it is blocked more than 30 seconds, it must be the worst case of threading
			except Empty as e:
				continue # released from the lock, restart a new workflow

			try:
				submitArgs = searchKeyword(task)
			except Exception as e:
				handleError(e, task)
				taskQueue.task_done() # search fails, it should be regared as a finished task

		if lastSubmitted: #check whether there is a submit before
			analysisResult = responseQueue.get() # if there is a submit, wait until there is an available analysis
		
		# submit the lastest sequence first before saving the analysis result
		# it saves the file writing time
		if submitArgs:
			try:
				submitToCRISPR(*submitArgs)
				newSubmit = True
			except Exception as e:
				handleError(e, task)
				newSubmit = False # submit fails
				taskQueue.task_done()
		else: # no submit request
			newSubmit = False

		if lastSubmitted:
			if analysisResult[1]: #make sure the response is available
				storeAnalysis(*analysisResult)

			taskQueue.task_done() # task isn't done until it is saved or no response

		lastSubmitted = newSubmit


def loadQueue(fileName):
	with open(fileName) as keywordsFile:
		i = 0
		for searchWord in keywordsFile:
			task = searchWord.rstrip()
			if task:
				taskQueue.put(task)
				i += 1

		global total
		total = i

		updateStatus()	

def workersStart(fileName, workerNum):

	print(progressTextTmp % 0, end = "")
	
	for i in range(workerNum):
		t = threading.Thread(target = worker)
		t.daemon = True
		t.start()

	loadQueue(fileName)
	taskQueue.join()

def writeHeaders():
	log.write(getTime() + " Gene Bot Started\n")
	analysisTableFile.write(rowTemplate.format( "Symbol",
												"Score",
												"Pool",
												"Sequence",
												"Direction",
												"Protein Bind"))
	log.flush()
	analysisTableFile.flush()

def doubleQuoteText(text):
	return '"' + text + '"'

def getTime():
	return time.strftime("%Y/%m/%d %H:%M:%S")

def makeOutputFolder():
	if not os.path.exists(outputDirName):
		print('see the output files in the %s folder' % doubleQuoteText(outputDirName))
		os.mkdir(outputDirName)

def readOptions():
	geneBotDes = """This is an internet robot.
					It searches and downloads the sequences of symbols specfied in the input file.
					Then, uploads these sequences to CRISPR server that analyzes them.
					Finally, downloads the analyses and parses them into a table as the output"""
	parser = argparse.ArgumentParser(description = geneBotDes, formatter_class = argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("filename", help = "specify a file to process")
	parser.add_argument("-o", "--range-offset", help = "specify the range offset of a sequence", type = int, default = 0)
	parser.add_argument("-l", "--seq-length", help = "specify the length of a sequence", type = int, default = 0)
	parser.add_argument("-e", "--email", help = "set the email when submit the sequence request to CRISPR server", default = "anonymous@whoareyou.com")
	parser.add_argument("-w", "--workers", help = "specify the thread(worker) number that process in parallel", type = int, default = 20)
	parser.add_argument("-i", "--watching-interval", help = "specify timer interval(in seconds) when watching the analysis result", type = int, default = 15)
	parser.add_argument("-t", "--watching-timeout", help = "specify the max waiting time(in minutes) when watching a analysis  ", type = int)
	parser.add_argument("-v", "--verbose", help = "output any process event when logging", action = "store_true")
	args = parser.parse_args()

	return (args.filename,
			args.range_offset,
			args.seq_length,
			args.email,
			args.workers,
			args.watching_interval,
			args.watching_timeout,
			args.verbose)


def main():
	global rangeOffset, seqLength, email, watchingInterval, watchingTimeout, logVerbose
	(fileName, rangeOffset, seqLength, email, workerNum, watchingInterval, watchingTimeout, logVerbose) = readOptions()

	print("Gene Bot Starting...")

	makeOutputFolder()
	if logVerbose:
		print("see live events in " + doubleQuoteText(logFilePath))

	global log, analysisTableFile
	with open(logFilePath, "w") as log, open(tableFilePath, "w") as analysisTableFile:
		writeHeaders()
		workersStart(fileName, workerNum)
		log.write(getTime() + " Process Complete. %d success(es), %d warning(s), and %d error(s)" % (completeCount, warningCount, errorCount))

	closeBackupFile()
	print("\nProcess Complete. see output table at " + doubleQuoteText(tableFilePath))
	print('%d warning(s), %d error(s), see %s for more details' % (warningCount, errorCount, doubleQuoteText(logFilePath)))


if __name__ == "__main__":
	main()
