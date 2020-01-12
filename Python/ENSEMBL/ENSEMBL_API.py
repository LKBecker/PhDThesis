#derived from ENSEMBLs example client, https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
#tested on Python 3.2.3 
#see also: http://rest.ensembl.org/documentation
import sys
import urllib.error
import urllib.parse
import urllib.request
import json
import time
import datetime

DEFAULT_SERVER = 'http://rest.ensembl.org'
DEFAULT_REQS_PER_SECOND = 13
DEFAULT_GENOMIC_REFERENCE = "GRCh37.p13"
MAX_SLICE_FOR_GENES = 500000
VERSION="1.3.0"

def ItemOrDefault(dict, key, default="[None found]"):
	if (key in dict.keys()): return dict[key] 
	else: return default

#Classes
class ENSEMBL_Object():
	def __init__(this, _dataDict):
		this.dataDict = _dataDict
				
	def __getattr__(self, item):
			if item in self.dataDict.keys(): return self.dataDict[item]
			else: 
				try:
					return object.__getattribute__(this, item)
				except:
					return "[None found]"

	def __str__(this):
		out = "\t".join([this.display_name, this.id, this.description, this.Parent, this.object_type, str(this.start), str(this.end), this.assembly_name])
		return out
		
class Transcript(ENSEMBL_Object):
	def __init__(this, transDict):
		this.dataDict = transDict
		
		if ((this.object_type!="Transcript") & (this.feature_type!="transcript")): raise ValueError("Attempted to create object of type 'transcript' from data of type '{}'".format(this.object_type))
			
		lookupStr = str('/lookup/id/' + this.Parent)
		this.parentData = client.perform_rest_action(lookupStr)
		if (this.parentData):
			this.parentData = str(this.parentData)
		else: this.parentData = ""
				
	def __str__(this):
			#Neigbor_Item_Name\tENSEMBL_ID\tNeighbor_Item Type\tBiological_Type\tSource\tNeighbor_Chromosome\tNeighbor_Start\tNeighbor_Stop\tNeighbor_Strand\tNeighbor_Assembly\Parent
			out = "\t".join([
				this.id, 
				this.transcript_id,
				this.feature_type,
				this.biotype,
				"N/A",
				this.source,
				str(this.seq_region_name),
				str(this.start),
				str(this.end),
				str(this.strand),
				this.assembly_name,
				this.Parent
			])
			return out
			
	def __getattr__(self, item):
		if item in self.dataDict.keys(): return self.dataDict[item]
		else: 
			try:
				return object.__getattribute__(this, item)
			except:
				return "[None found]"

class Gene(ENSEMBL_Object):
	def __init__(this, dataDict):
		this.dataDict = dataDict
		
		# try: assert 'object_type' in this.dataDict
		# except AssertionError: this.dataDict['object_type'] = "[None found]"
		
		if ((this.object_type != "Gene") & (this.feature_type != "gene")): raise ValueError("Attempted to create object of type 'gene' from data of type '{}'".format(this.object_type))
	
		try: assert  'description' in this.dataDict
		except AssertionError: this.dataDict['description'] = "[None found]"
		if(this.dataDict['description']==None): this.dataDict['description']="[None found]"

	def __getattr__(self, item):
		if item in self.dataDict.keys(): return self.dataDict[item]
		else: 
			try:
				return object.__getattribute__(this, item)
			except:
				return "[None found]"

	def __str__(this):
			#Neigbor_Item_Name\tENSEMBL_ID\tNeighbor_Item Type\tBiological_Type\tSource\tNeighbor_Chromosome\tNeighbor_Start\tNeighbor_Stop\tNeighbor_Strand\tNeighbor_Assembly\Parent
			out = "\t".join([
				this.display_name, 
				this.gene_id, 
				this.feature_type, 
				this.biotype, 
				this.description, 
				this.source, 
				str(this.seq_region_name), 
				str(this.start), 
				str(this.end), 
				str(this.strand),
				this.assembly_name,
				"(is Parent)"
			])
			return out

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
			req = urllib.request.Request(this.server + endpoint, data=json_data, headers = {'content-type':'application/json'})
			response = urllib.request.urlopen(req)
			content = response.read().decode('utf-8')
			if content: JSONdata = json.loads(content)	#should always be JSON
			this.req_count += 1

		except urllib.error.HTTPError as e: # check if we are being rate limited by the server
			if e.code == 429:
				if 'Retry-After' in e.headers:
					retry = e.headers['Retry-After']
					time.sleep(float(retry))
					this.perform_rest_action(endpoint, params)	#recurse this after specified retry period
			else:
				sys.stderr.write('Request failed for {0}: {1.code} {1.reason}\n'.format(endpoint, e))
		return JSONdata

	def get_variants(this, species, symbol):
		genes = this.perform_rest_action(
			'/xrefs/symbol/{0}/{1}'.format(species, symbol), 
			params={'object_type': 'gene'}
		)
		if genes:
			stable_id = genes[0]['id']
			variants = this.perform_rest_action(
				'/overlap/id/{0}'.format(stable_id),
				params={'feature': 'variation'}
			)
			return variants
		return None
		
	def get_SNP(this, species, symbol):
		SNP = this.perform_rest_action('/variation/{0}/{1}'.format(species, symbol))
		if SNP:
			DNAData = SNP['mappings'][0]
			data = []
			location = DNAData['location'].split(':')
			location[1] = location[1].split("-")[0]
			data.append(SNP["name"])
			data.append(DNAData['allele_string'])
			data.append(SNP['ancestral_allele'])
			data.append(SNP['minor_allele'])
			data.append(DNAData['assembly_name'])
			data.append(DNAData['strand'])
			data.append(location[0])
			data.append(location[1])
			if(len(SNP["synonyms"])>0): data.append("\t".join(SNP["synonyms"]))
			else: data.append("[None found]")
			return data
		return [symbol, "[None found]", "[None found]", "[None found]", "[None found]", "[None found]", "[None found]", "[None found]", "[None found]"]
		
	def get_Feature(this, symbol):
		data = this.perform_rest_action('/lookup/id/{0}'.format(symbol))
		if data: 
			if symbol[:4]=="ENSG":
				_data = Gene(data)
			elif symbol[:4]=="ENST":
				_data = Transcript(data)
			elif symbol[:4]=="ENSO": #REDO - ENSOARG doesn't mean ENS OBJ but ENS Ovis ARies Gene
				_data = ENSEMBL_Object(data)				
			else:
				return "Cannot parse unknown type " + symbol[:4]
			return _data
		return None
			
	#/overlap/region/human/7:140424943-140624564?feature=gene;feature=transcript;feature=cds;feature=exon;content-type=application/json
	def get_Features_around_SNP(this, species, chromosome, start, stop, features):
		if stop < start:
			a = start
			start = stop
			stop = a
			del a
		
		#todo - redo??
		if (abs(stop-start) > MAX_SLICE_FOR_GENES): 
			print("Range too high! Adjusting from both sides...")
			a = (abs(stop-start)-MAX_SLICE_FOR_GENES)/2
			start = int(start + a)
			stop = int(stop - a)
			
		#query = '/overlap/region/{}/{}:{}-{}?feature=regulatory'.format(species, chromosome, start, stop)
		#Enum(band, gene, transcript, cds, exon, repeat, simple, misc, variation, somatic_variation, structural_variation, somatic_structural_variation, constrained, regulatory, segmentation, motif, chipseq, array_probe)
		featurestring = ";".join(map(lambda x: "feature="+x, features))
		query = '/overlap/region/{}/{}:{}-{}?{}'.format(species, chromosome, start, stop, featurestring)
		#print(query)
		Features = this.perform_rest_action(query)
		if Features:
			data = []
			for feature in Features:
				#print(feature)
				if(feature['feature_type']=="transcript"): 
					_feature = Transcript(feature)
				elif (feature['feature_type']=="gene"): _feature = Gene(feature)
				else:
					print("WARNING: Unknown feature type:" + feature['feature_type'])
					continue
				data.append(_feature)
			return data
		return None
		
		#/overlap/region/human/7:140424943-140624564?feature=gene;feature=transcript;feature=cds;feature=exon;content-type=application/json
	def get_Overlapping_Features(this, species, chromosome, start, stop, features):
		if stop < start:
			a = start
			start = stop
			stop = a
			del a
					
		#query = '/overlap/region/{}/{}:{}-{}?feature=regulatory'.format(species, chromosome, start, stop)
		#Enum(band, gene, transcript, cds, exon, repeat, simple, misc, variation, somatic_variation, structural_variation, somatic_structural_variation, constrained, regulatory, segmentation, motif, chipseq, array_probe)
		query = '/overlap/region/{}/{}:{}-{}'.format(species, chromosome, start, stop)
		query += "?"
		for feature in features: query += "feature=" + feature + ";"
		#print(query)
		Features = this.perform_rest_action(query)
		if Features:
			data = []
			for feature in Features:
				#print(feature)
				if(feature['feature_type']=="transcript"): 
					_feature = Transcript(feature)
				elif (feature['feature_type']=="gene"): _feature = Gene(feature)
				else:
					print("WARNING: Unknown feature type:" + feature['feature_type'])
					continue
				data.append(_feature)
			return data
		return None	
		
	def remap_to_other_Assembly(this, species, oldAssembly, newAssembly, oldChromosome, oldStart, oldEnd, oldStrand = 1):
		endpoint = '/map/{}/{}/{}:{}..{}:{}/{}'.format(species, oldAssembly, oldChromosome, oldStart, oldEnd, oldStrand, newAssembly)
		remap = this.perform_rest_action(endpoint, {'coord_system':'chromosome', 'target_coord_system':'chromosome'})['mappings']
		if(len(remap)>=1):
			remap = remap[0]['mapped']
			#new assembly	chromosome	start	stop	strand
			result = "\t".join(map(lambda y: str(y), [remap["assembly"], remap["seq_region_name"],remap["start"],remap["end"],remap["strand"]]))
			return result
		return None
		
	#def get_Gene_location(this, geneID):
	#	endpoint = '/lookup/id/{}'.format(geneID)
	#	geneLoc = this.perform_rest_action(endpoint)
	#	if(len(geneLoc)>=1):
	#		#TODO PROCESS OUTPUT
	#		#return geneLocProcessed
	#	return None
		
	# /vep/human/id/COSM476?content-type=application/json
	def get_SNP_Consequence(this, species, rsID):
		endpoint='/vep/{}/id/{}'.format(species, rsID.strip())
		VEPdata = this.perform_rest_action(endpoint)
		if(VEPdata):
			return VEPdata[0]
		return None
		
	def get_condensed_human_homolog(this, geneID):
		#;target_species=human;format=condensed
		endpoint='/homology/id/{}'.format(geneID)
		HomologData = this.perform_rest_action(endpoint, params={'target_species':'human', 'format':'condensed'})['data']
		if (len(HomologData) > 0): HomologData = HomologData[0]['homologies']
		else: return("None Found")
		if(len(HomologData) > 0):
			result = HomologData[0]['id']
			return(result)
		else: return('None Found')
		
	#currently unused #TODO
	def get_sequence_from_region(this, chromosome, start, stop, species="human", strand="1", mask_feature=0):
		if((stop-start)>10e+6): stop = start + 10e+6
		endpoint='/sequence/region/{}/{}:{}..{}:{}'.format(species, chromosome, start, stop, strand)
		params={'content-type':'text/plain', 'mask_feature':mask_feature}
		SequenceData = this.perform_rest_action(endpoint, params)
		if (SequenceData): return SequenceData
		else: return None
		
	def get_sequence_from_identifier(this, ID, seqType):
		if seqType not in ["genomic", "cds", "cdna", "protein"]: raise ArgumentException("Error in get_sequence_from_identifier(): type must be one of [genomic/cdna/cda/protein], is " + seqType)
		endpoint='/sequence/id/{}'.format(ID)
		params={'type':seqType}
		if seqType == "protein" : params["multiple_sequences"]=1
		SequenceData = this.perform_rest_action(endpoint, params=params)
		if(SequenceData): return SequenceData
		else: return None
		
	def get_r2_value(this, snp1, snp2, species, population=None):
		if snp1=='' or snp2=='': return None
		endpoint='/ld/{}/pairwise/{}/{}'.format(species, snp1, snp2)
		params={}
		if population is not None: params['population']=population
		R2Data = this.perform_rest_action(endpoint, params=params)
		if(R2Data): return R2Data
		else: 
			return None
	
#Functions
def batch_SNPs(species="human"):
	#expected input: single column, rsID according to dbSNP
	# client = EnsemblRestClient()
	with open("./input SNP.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL SNP Data.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: SNP data retrieval\nParameters:\nSpecies\t{}\n\n".format(VERSION, runtime,species))
			out = "SNPID\tAlleleString\tAncestralAllele\tMinorAllele\tAssembly\tStrand\tChromosome\tPosition\tSynonyms\n"
			print(out)
			output.write(out)
			for line in inputfile:
				data = client.get_SNP(species, line.strip())
				data = map(lambda x: str(x), data)
				out = "\t".join(data)
				print(out)
				output.write(out+"\n")
	return
	
def batch_SNP_Consequences(species="human"):
	# client = EnsemblRestClient()
	with open("./input SNP.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL VEP output.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: Variant-Effect Prediction query\nParameters:\nSpecies\t{}\n\n".format(VERSION, runtime,species))
			out = "\nSNPID\tAssociatedGeneID\tAssociatedTranscriptID\tAssociatedGeneName\tAssociatedGeneType\tImpactRating\tVariantAllele\tConsequenceTerms\n"
			print(out)
			output.write(out)
			for line in inputfile:
				data = client.get_SNP_Consequence(species, line)
				SNP = data["id"]
				MSC = data["most_severe_consequence"]
				if ("transcript_consequences" in data.keys()):
					for AssocDataDict in data["transcript_consequences"]:
						GeneID 			= ItemOrDefault(AssocDataDict, "gene_id")
						TranscriptID	= ItemOrDefault(AssocDataDict, "transcript_id")
						GeneName 		= ItemOrDefault(AssocDataDict, "gene_symbol")
						GeneType		= ItemOrDefault(AssocDataDict, "biotype")
						Impact 			= ItemOrDefault(AssocDataDict, "impact")
						VariantAllele 	= ItemOrDefault(AssocDataDict, "variant_allele")
						Consequences  	= "/".join(ItemOrDefault(AssocDataDict, "consequence_terms", ["[None found]"]))
						out="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(SNP, GeneID, TranscriptID, GeneName, GeneType, Impact, VariantAllele, Consequences)
						print(out)
						output.write(out+"\n")
				if("intergenic_consequences" in data.keys()):
					for AssocDataDict in data["intergenic_consequences"]:
						GeneID 			= "(Intergenic; no Gene)"
						TranscriptID	= "(Intergenic; no Transcript)"
						GeneName 		= "(Intergenic; no Gene)"
						GeneType		= "(Intergenic; not expressed)"
						Impact 			= ItemOrDefault(AssocDataDict, "impact")
						VariantAllele 	= ItemOrDefault(AssocDataDict, "variant_allele")
						Consequences  	= "/".join(ItemOrDefault(AssocDataDict, "consequence_terms", ["[None found]"]))
						out="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(SNP, GeneID, TranscriptID, GeneName, GeneType, Impact, VariantAllele, Consequences)
						print(out)
						output.write(out+"\n")

def batch_SNP_Contextualise(species="human", features=["gene"]):
#expected input SNPID	Chromosome	Location (single number)
	with open("./input feature search.txt" , 'r') as input:
		with open("./{} EMSEMBL Batch Feature Search output.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: In-Range search\nParameters:\nSearch range\t{}bp either side\nFeature type(s)\t{}\nSpecies\t{}\n\n".format(VERSION, runtime, MAX_SLICE_FOR_GENES/2, ";".join(features), species))
			out = "Base_SNP\tSNP_Chromosome\tSNP_Position\tNeigbor_Item_Name\tENSEMBL_ID\tNeighbor_Item Type\tBiological_Type\tSource\tChromosome\tStart\tStop\tAssembly\n"
			print(out)
			output.write(out)
			for line in input:
				line = line.split("\t")
				baseSNP = 		line[0]
				chromosome =	int(line[1])
				start = 		int(max(1, int(line[2]) - (MAX_SLICE_FOR_GENES/2)))
				stop = 			int(int(line[2]) + (MAX_SLICE_FOR_GENES/2))
				data = 			client.get_Features_around_SNP(species, chromosome, start, stop, features)
				#print(data)
				if data is not None:
					for item in data:
						out="\t".join(list(map(lambda x: str(x), [baseSNP, chromosome, line[2].rstrip(), item])))
						print(out)
						output.write(out+"\n")
	return
	
def batch_Feature_Overlap(species="human", features=["gene"]):
	#expected input RefID	Chromosome	Start Stop
	with open("./input feature search.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL Batch Feature Search output.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: In-Range search\nParameters:\nSearch range\t{}bp either side\nFeature type(s)\t{}\nSpecies\t{}\n\n".format(VERSION, runtime, MAX_SLICE_FOR_GENES/2, ";".join(features), species))
			out = "Base_Feature\tBase_Chromosome\tBase_Start\tBase_Stop\tNeigbor_Item_Name\tENSEMBL_ID\tNeighbor_Item_Type\tBiological_Type\tName\tSource\tNeighbor_Chromosome\tNeighbor_Start\tNeighbor_Stop\tNeighbor_Strand\tNeighbor_Assembly\tParent\n"
			print(out)
			output.write(out)
			for line in inputfile:
				line = line.split("\t")
				baseFeature = 		line[0]
				chromosome =	str(line[1])
				start = 		int(line[2])
				if len(line) < 4: 
					stop = 		start + 1
				else:
					stop = 		int(line[3])
				search_start = 	int(max(1, (start - (MAX_SLICE_FOR_GENES/2))))
				search_stop = 	int(stop + (MAX_SLICE_FOR_GENES/2))
				data = 			client.get_Overlapping_Features(species, chromosome, search_start, search_stop, features)
				#print(data)
				if data is not None:
					for item in data:
						out="\t".join(list(map(lambda x: str(x), [baseFeature, chromosome, start, stop, item])))
						print(out)
						output.write(out+"\n")
	return
	
def batch_Characterise():
	#expected input: ENSEMBL_ID
	with open("./input characterise.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL Batch Feature contextualisation output.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: Characterise\n\n".format(VERSION, runtime, ))
			out = "Feature_ID\tName\tAvailable Data\n"
			print(out)
			output.write(out)
			lineCount = 0
			for line in inputfile:
				lineCount+=1
				line = line.split("\t")
				ID = 			line[0].strip()
				data = 			client.get_Feature(ID)
				if data is not None:
					out="\t".join(list(map(lambda x: str(x), [ID, data])))
					print(out)
					output.write(out+"\n")
					#if(lineCount%100==0): print("Processing line {}".format(lineCount))
	return
	
def batch_remap(species="human"):
	#expected input: tab-separated, input assembly (GrCH[X]), chromosome, start, stop, output assembly (GrCh[X])
	# client = EnsemblRestClient()
	with open("./input remap.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL Batch Remap output.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: batch remap.\nParameters:\nSpecies\t{}\n\n".format(VERSION, runtime, species))
			out = "Input Assembly\tChromosome\tStart\tEnd\tStrand\tOutput Assembly\tChromosome\tStart\tEnd\tStrand\n"
			print(out)
			output.write(out)
			for line in inputfile:
				line = line.split("\t")
				InputAssembly =		line[0]
				InputChr =			line[1]
				InputStart = 		line[2]
				InputStop =			line[3]
				RequestedAssembly =	line[4]
				remap = client.remap_to_other_Assembly(species, InputAssembly, RequestedAssembly, InputChr, InputStart, InputStop)	
				out="\t".join(list(map(lambda x: str(x), [InputAssembly, InputChr, InputStart, InputStop, "1", remap])))
				print(out)
				output.write(out+"\n")	

def batch_human_homolog():
	#expected input: tab-separated, input assembly (GrCH[X]), chromosome, start, stop, output assembly (GrCh[X])
	# client = EnsemblRestClient()
	with open("./input human homolog.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL Human Homolog.txt".format(runtime), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: batch human homolog\n\n".format(VERSION, runtime))
			out = "Input ID\tHuman homolog\n"
			output.write(out)
			print(out)
			for line in inputfile:
				ID = line.strip()
				remap = client.get_condensed_human_homolog(ID)	
				out="{}\t{}".format(ID, remap)
				print(out)
				output.write(out+"\n")	

def batch_get_sequence(type="gene"):
	#expected input: ENSEMBL IDs, one ID per line. ID determines species.
	# client = EnsemblRestClient()
	with open("./input get sequence.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL {} Sequence Retrieval.txt".format(runtime, type), 'w') as output:
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: batch sequence retrieval\n\n".format(VERSION, runtime))
			out = "Input ID\tProteinID\tSequence\n"
			output.write(out)
			print(out)
			for line in inputfile:
				ID = line.strip()
				seqData = client.get_sequence_from_identifier(ID, type)
			# for subseq in SequenceData:
				# if "error" in subseq.keys(): return None
				# if subseq["molecule"] != seqType: 
					# print("got: "+ subseq["molecule"] + ", requested: " + seqType)
					# return None
			# else: return SequenceData
				if seqData:	
					for seq in seqData:
						if (seq["molecule"] == type):
							out="{}\t{}\t{}".format(ID, seq["id"], seq["seq"])
							
						else: 
							out="{}\t{}\t{}".format(ID, seq["id"], "Warning: Mismatch. " + type + " requested, got " + seq["molecule"] + " instead.")
				else:
					out="{}\t0\tNone Found".format(ID)
				print(out)
				output.write(out+"\n")
						
def batch_SNP_r2_retrieval(species="human", population="1000GENOMES:phase_3:GBR"):
	with open("./input rsquared retrieval.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL r2 value retrieval.txt".format(runtime), 'w') as output:
			out = ""
			output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: batch SNP-SNP r2 value retrieval\nParameters:\nspecies\t{}\npopulation\t{}\n\n".format(VERSION, runtime, species, population))
			output.write(out)
			print(out)
			line1 = inputfile.readline()
			line2 = inputfile.readline()
			if not inputfile.readline() == '': raise Exception("File appears to have wrong format")
			SNPs1 = line1.split(" ")
			SNPs1 = list(filter(lambda y: y!='', map(lambda x: x.strip(), SNPs1)))
			SNPs2 = line2.split(" ")
			SNPs2 = list(filter(lambda y: y!='', map(lambda x: x.strip(), SNPs2)))
			PairsDone=[]
			output.write("\t" + "\t".join(SNPs2) + "\n")
			for i, snp in enumerate(SNPs1):
				output.write(snp)
				for _snp in SNPs2:
					if (((snp+_snp) in PairsDone) or ((_snp+snp) in PairsDone)): 
						output.write("\t")
						continue						
					else: 
						PairsDone.append((snp+_snp))
						r2 = client.get_r2_value(snp, _snp, species, population)
						if r2 is not None:
							r2 = r2[0]['r2']
							output.write("\t" + r2)
							print(snp, "*", _snp, ":", r2)
						else:
							output.write("\t")
							print(snp, "*", _snp, ": NA")
					output.write("\n")
					print("Progress: {}%".format((i/len(SNPs1))*100))
				
def batch_SNP_MAF_retrieval(species="human"):
	with open("./input Get MAF.txt" , 'r') as inputfile:
		SNPs = inputfile.read().split("\n")
	SNPs = list(filter(lambda x: len(x) > 0, SNPs))
	endpoint='/variation/{}'.format(species)
	MAFDict = {'ids':SNPs, 'pops':1}
	MAFData = client.perform_rest_action(endpoint, data=MAFDict)
	with open("./{} EMSEMBL MAF retrieval.txt".format(runtime), 'w') as output:
		output.write("Running ENSEMBL API client, version {}.\nStartup time: {}\nSubroutine: batch SNP MAF retrieval\nParameters:\nSpecies\t{}\n\n".format(VERSION, runtime, species))
		out = "rsID\tMinor Allele\tMAF\tAncestral allele\n"
		output.write(out)
		print(out)	
		for key in MAFData.keys():
			out = "{}\t{}\t{}\t{}\n".format(key, str(MAFData[key]['minor_allele']), str(MAFData[key]['MAF']), str(MAFData[key]['ancestral_allele']))
			print(out)
			output.write(out)

def batch_UpstreamSequences(species="sheep", upstreamBp=5000):
	with open("./input get Upstream 5k.txt" , 'r') as inputfile:
		IDs = inputfile.read().split("\n")
		with open("./{} EMSEMBL upstream sequence retrieval.txt".format(runtime), 'w') as outputCorr:
			with open("./{} EMSEMBL upstream sequence retrieval errors.txt".format(runtime), 'w') as outputErr:
				for ID in IDs:
					data = client.get_Feature(ID.strip())
					if type(data) != str:
						try: Coord_Chr	= data.seq_region_name
						except:
							errStr="Cannot locate chromosome information for gene {}. Please try manually.".format(ID)
							print(errStr)
							outputErr.write(errStr)
							continue
						
						if ((data.start - (upstreamBp +1)) < 1):
							errStr = "Gene '{}' appears to not have {}bp of upstream space. Check if region '{}' is a satellite or try manually.\n".format(ID, upstreamBp, Coord_Chr)
							print(errStr)
							outputErr.write(errStr)
							continue
							
						Coord_Start = data.start - (upstreamBp +1)
						Coord_End	= data.start
						try: Coord_Assem = data.assembly_name
						except: Coord_Assem = "ERROR"
						print("Found entry for ID '{}'. Attempting to retrieve {}:{}-{} ({})...".format(ID, Coord_Chr, Coord_Start, Coord_End, Coord_Assem))
						Seq_Data = client.get_sequence_from_region(Coord_Chr, Coord_Start, Coord_End, species) #strand="1"
						
						if Seq_Data and 'seq' in Seq_Data: outputCorr.write("{}\t{}\n".format(ID, Seq_Data['seq']))
						else: 
							errStr = "Failed to retrieve sequence data for ID '{}', location {}:{}-{}.\n".format(ID, Coord_Chr, Coord_Start, Coord_End)
							print(errStr)
							outputErr.write(errStr)
							
					else: 
						errStr = "Failed to retrieve data for ID '" + ID + "'.\n"
						print(errStr)
						outputErr.write(errStr)
			print(" ---- Complete.")

def batch_get_sequence_by_coords(species="human", mask=0, strand=1):
	with open("./input get seq from coords.txt" , 'r') as inputfile:
		with open("./{} EMSEMBL sequence retrieval.txt".format(runtime), 'w') as outputCorr:
			Seqs = inputfile.read().split("\n")
			for Seq in Seqs:
				name, Coord_Chr, Coord_Start, Coord_End = Seq.split("\t")
				print("Attempting to retrieve sequence for '{}' {}:{}-{} ({})...".format(name, Coord_Chr, Coord_Start, Coord_End, strand))
				Seq_Data = client.get_sequence_from_region(Coord_Chr, Coord_Start, Coord_End, species)
				outputCorr.write("")
						
if __name__ == "__main__":
	global client
	global runtime 
	runtime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
	client = EnsemblRestClient()
	#batch_SNPs()
	#batch_SNP_Consequences()
	#batch_Feature_Overlap("human", ["gene"])
	#batch_Feature_Overlap("sheep", ["cds"])
	#batch_SNP_Contextualise()
	#batch_remap("human")
	#batch_Characterise()
	#batch_human_homolog()
	#batch_get_sequence("protein")
	#batch_SNP_r2_retrieval()
	#batch_SNP_MAF_retrieval()
	#batch_UpstreamSequences()
	batch_get_sequence_by_coords("sheep")
