from os import listdir
from os.path import isfile, join
import re
import datetime

HiSat2ResultFile = re.compile("(?P<reads>\d+) reads; of these:\s+(?P<nPairedReads>\d+) \(\d{1,3}\.\d{2}%\) were paired; of these:\s+(?P<PRAlignedConcordant0Times>\d+) \(\d{1,3}\.\d{2}%\) aligned concordantly 0 times\s+(?P<AlignedConcord1Time>\d+) \(\d{1,3}\.\d{2}%\) aligned concordantly exactly 1 time\s+(?P<PRConcordManyTime>\d+) \(\d{1,3}\.\d{2}%\) aligned concordantly >1 times\s+----\s+\d+ pairs aligned concordantly 0 times; of these:\s+(?P<PRConcord0Discord1>\d+) \(\d{1,3}\.\d{2}%\) aligned discordantly 1 time\s+----\s+(?P<PRNeverAligned>\d+) pairs aligned 0 times concordantly or discordantly; of these:\s+\d+ mates make up the pairs; of these:\s+(?P<NAPNAMates>\d+) \(\d{1,3}\.\d{2}%\) aligned 0 times\s+(?P<NAP1AlignMates>\d+) \(\d{1,3}\.\d{2}%\) aligned exactly 1 time\s+(?P<NAPManyAlignMates>\d+) \(\d{1,3}\.\d{2}%\) aligned >1 times\s+(?P<overallAlignmentRate>\d{1,3}\.\d{2})% overall alignment rate")

if __name__ == "__main__":
	runtime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
	Files = listdir()
	Files = filter(lambda x: x[-3:] not in [".py", "pyw", "txt"], Files)
	Files = list(filter(lambda x: "." in x, Files))
	
	Results = []
	with open("./{} HiSat2 Run data.txt".format(runtime), 'w') as output:
		output.write("ID\tReads\tOfWhichPaired\tAlignedConcord0Times\tAlignedConcord1Time\tAlignedConcordManyTime\tConcord0Discord1\tNeverAlignPairs\tNAP.NAMates\tNAP.1AlignMates\tNAP.ManyAlignMates\tOverallAlignmentRate\n")
		for file in Files:
			ID = int(file.split(".")[-1])
			with open(file, 'r') as io:
				FileData = HiSat2ResultFile.match(io.read())
			FileData = list(map(lambda x: float(x), FileData.groups()))
			output.write("{}\t{}\n".format(ID, "\t".join(map(lambda y: str(y), FileData))))
			