from urllib.request import urlopen
from queue import Queue
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

crisprHost = "http://crispr.mit.edu"
gbDownloadUrl = "/export/guides_gb/{0}"

tableDataRegex = re.compile(r'protein_bind\s+(.+?)\n\s+/bound_moiety="CRISPR\s([\w\-]+?)\sguide\s\#\d+\s+([A-Z]+)".+?"score":\s"(\d+%)"', re.DOTALL)


outputDirName = "output"

def makePath(fileName):
	return outputDirName + "/" + fileName


tableFilePath = makePath("recovered_analysis_table.txt")

total = -1
completeCount = 0

taskQueue = Queue()

progressTextTmp = "\rrecovered %d analyses"
refreshSpaceLen = 0

analysisTableFile = None
tableFileLock = threading.Lock()
rowTemplate =  "{}" + "\t{}" * 5 + "\n"

watchingInterval = 0

def storeTable(table):
	with tableFileLock:
		for row in table:
			analysisTableFile.write(rowTemplate.format( row["symbol"],
														row["score"],
														row["pool"],
														row["sequence"],
														row["direction"],
														row["protein_bind"]))

		analysisTableFile.flush()

def parseAnalysis(symbol, analysis):
	table = []
	for match in tableDataRegex.finditer(analysis):
		row = { "symbol": symbol,
		 		"protein_bind": match.group(1),
		 		"direction": match.group(2),
		 		"sequence": match.group(3),
		 		"score": match.group(4),
		 		"pool": "NUCLEAR"}

		table.append(row)

	return table

def storeAnalysis(symbol, response):
	
	table = parseAnalysis(symbol, getContentAndClose(response))
	table.sort(key = operator.itemgetter("score"), reverse = True)
	storeTable(table)

	finishTask(symbol)

def watchAnalysis(symbol, jobKey):
	resultUrl = crisprHost + gbDownloadUrl.format(jobKey)
	statusCode = 500
	response = None
	global analysisText
	while statusCode == 500:
		try:
			response = urlopen(resultUrl)
			statusCode = 200 #sever process complete
			break

		except Exception:
			updateStatus("still analyzing symbol " + doubleQuoteText(symbol))
			pass # force waiting until the analysis is fetched


				 
		time.sleep(watchingInterval)

	storeAnalysis(symbol, response)
	

def getContentAndClose(response, amt = None):
	content = response.read(amt).decode()
	response.close()
	return content

def finishTask(symbol):
	global completeCount
	completeCount = completeCount + 1
	updateStatus("symbol " + doubleQuoteText(symbol) + " is stored")


def updateStatus(notification = ""):
	progressText = progressTextTmp % completeCount
	percentageText = ""

	if total > 0:
		percentageText = " (%.2f%%)" % (completeCount * 100 / total)

	global refreshSpaceLen
	progressText += percentageText + ((": " + notification) if notification else "")
	stdout.write("\r" + " " * refreshSpaceLen)
	stdout.write(progressText)
	stdout.flush()
	refreshSpaceLen = len(progressText)

def worker():
	while True:			
		task = taskQueue.get()
		watchAnalysis(*task)
		taskQueue.task_done()


def loadQueue(fileName):
	with open(fileName) as backupFile:
		i = 0
		for item in backupFile:
			task = item.rstrip().split()
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
	analysisTableFile.write(rowTemplate.format( "Symbol",
												"Score",
												"Pool",
												"Sequence",
												"Direction",
												"Protein Bind"))
	analysisTableFile.flush()

def doubleQuoteText(text):
	return '"' + text + '"'

def makeOutputFolder():
	if not os.path.exists(outputDirName):
		print('see the output files in the %s folder' % doubleQuoteText(outputDirName))
		os.mkdir(outputDirName)

def readOptions():
	des = """This is the recovery program that read backup file and fetch the finished analysis automatically.
			Unlike Gene Bot, it ignores errors and try again when fetching the analysis."""
	parser = argparse.ArgumentParser(description = des, formatter_class = argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("filename", help = "specify a file to process")
	
	parser.add_argument("-w", "--workers", help = "specify the thread(worker) number that process in parallel", type = int, default = 20)
	parser.add_argument("-i", "--watching-interval", help = "specify timer interval (in seconds) when watching the analysis result", type = int, default = 15)
	
	args = parser.parse_args()

	return (args.filename,
			args.workers,
			args.watching_interval)


def main():
	global watchingInterval, analysisTableFile
	(fileName, workerNum, watchingInterval) = readOptions()

	print("Recovery Starting...")
	makeOutputFolder()

	with open(tableFilePath, "w") as analysisTableFile:
		writeHeaders()
		workersStart(fileName, workerNum)

	print("\nProcess Complete. see output table at " + doubleQuoteText(tableFilePath))

if __name__ == "__main__":
	main()
