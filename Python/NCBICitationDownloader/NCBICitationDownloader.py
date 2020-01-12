#derived from ENSEMBLs example client, https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
#tested on Python 3.2.3 
#see also: https://api.ncbi.nlm.nih.gov/lit/ctxp
import re
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import time
import datetime

#/lit/ctxp/v1/pmc/?format=ris&id=1301148
#/lit/ctxp/v1/pubmed/?id=10024300
DEFAULT_SERVER = 'https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/'
DEFAULT_REQS_PER_SECOND = 5
VERSION="0.0.1"
PROJECT="NCBILoader"
EMAIL="a@a.com"

NCBIUrl = re.compile(r'^(http(s)?://)?(www\.)?(ncbi\.nlm\.nih\.gov/)?')
NCBIPaper = re.compile(r'^(?P<database>pubmed|pmc)/(articles/)?(?P<id>(PMC)?\d{6,})')

#Contains actual processing!
class EnsemblRestClient():
	def __init__(this, server=DEFAULT_SERVER, reqs_per_sec=DEFAULT_REQS_PER_SECOND):
		this.server = server
		this.reqs_per_sec = reqs_per_sec
		this.req_count = 0
		this.last_req = 0

	def perform_rest_action(this, endpoint, params=None, data=None):
		if params is not None: endpoint += '?' + urllib.parse.urlencode(params)
		JSONdata = None
		
		# check if we need to rate limit ourselves
		if this.req_count >= this.reqs_per_sec:
			delta = time.time() - this.last_req
			if delta < 1:
				time.sleep(1 - delta)			#sleep until limit normal again
			this.last_req = time.time()
			this.req_count = 0					#reset limit

		try:
			if (data is not None): json_data = json.dumps(data).encode('utf8')
			else: json_data = None
			req = urllib.request.Request(this.server + endpoint, data=json_data, headers = {'content-type':'application/json', 'User-Agent':PROJECT})
			response = urllib.request.urlopen(req)
			content = response.read().decode('utf-8', 'replace')
			this.req_count += 1

		except urllib.error.HTTPError as e: # check if we are being rate limited by the server
			if e.code == 429:
				if 'Retry-After' in e.headers:
					retry = e.headers['Retry-After']
					time.sleep(float(retry))
					this.perform_rest_action(endpoint, params)	#recurse this after specified retry period
			else:
				sys.stderr.write('Request failed for {0}: {1.code} {1.reason}\n'.format(endpoint, e))
				return None
		return content #JSONdata

	def get_citation(this, database, id, citationFormat):
		citation = this.perform_rest_action( '/{0}'.format(database), params={'format': citationFormat, 'id': id} )
		return citation
		
	
#Functions
def batch_get_citations_as_RIS(citationFormat="ris"):
	#expected input: single column, rsID according to dbSNP
	
	outputFolder = "./" + runtime + "_CitationDownload"
	try:
		os.mkdir(outputFolder)
	except:
		print("Could not create folder for citation files. Please ensure you have permissions and retry.")
		return

	with open("./input_citations.txt" , 'r') as inputfile:
		for line in inputfile:
			if (NCBIUrl.match(line) is None): 
				print("Invalid line: " + line)
				continue
			line = 	NCBIUrl.sub("", line)
			Match = re.match(NCBIPaper, line)
			if Match:
				db = Match.group("database")
				id = Match.group("id")
				if (db.lower() == "pmc"): id = id[3:] #cut off the "PMC" part of the ID
				
				data = client.get_citation(db, id, citationFormat)
				if data is not None:
					FirstAuthor = re.search("^AU\s+- (.+),", data, flags=re.M).group(1)
					Year = re.search("^Y1  - (\d{4})", data, flags=re.M).group(1)
					FileName="{0}{1}{2}{3}".format(outputFolder, os.path.sep, Year, FirstAuthor)
					
					if (os.path.isfile("{0}.{1}".format(FileName, citationFormat))):
						Counter = 1
						while os.path.isfile("{0}({1}).{2}".format(FileName, Counter, citationFormat)): Counter += 1
						FileName = "{0}({1}).{2}".format(FileName, Counter, citationFormat)
					else: FileName = "{0}.{1}".format(FileName, citationFormat)
					#print(FileName)
					
					with open(FileName, 'w', encoding="utf-8") as output: output.write(data)
					print("Citation for {0} ({1}) saved to file.".format(FirstAuthor, Year))
					
				else:
					print("Error attempting to retrieve citation for {0} {1}".format(db, id))
					
			else:
				print("Error attempting to extract database and ID from URL " + line)
	return
	
if __name__ == "__main__":
	global client
	global runtime 
	runtime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
	client = EnsemblRestClient()
	batch_get_citations_as_RIS()
