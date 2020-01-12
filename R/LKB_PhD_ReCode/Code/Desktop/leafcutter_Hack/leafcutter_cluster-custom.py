#Exmaination of the BAM files of aligned reads has shown that introns, despite support from more than -m 10 reads, disappear in the clustering step
#This implies some kind of overactive cleaning mechanism
#Thus, we're mucking about with the clustering file. GitHub, GO!

#!/usr/bin/env python
import sys
import tempfile
import math
import os
import gzip
import shutil
import datetime

BIN1 = ['2', 117267199, 117314475]

def merge_junctions(): 
	''' function to merge junctions '''
	fnameout = "%s/%s"%(rundir,outPrefix)
	flist = "%s/%s_sortedlibs"%(rundir, outPrefix)
	
	lsts = []
	for ln in open(flist):
		lsts.append(ln.strip())
	print("merging %d junction files...\n"%(len(lsts)))
	
	# Change 300 if max open file is < 300
	N = min([300, max([100, int(len(lsts)**(0.5))])])

	tmpfiles = []
	while len(lsts) > 1: 
		clst = []
		
		for i in range(0,(len(lsts)/N)+1): 
			lst = lsts[N*i:N*(i+1)]
			if len(lst) > 0:
				clst.append(lst)
		lsts = []
	
		for lst in clst:
			if len(lst) == 0: continue
			tmpfile = tempfile.mktemp()
			os.mkdir(tmpfile)
			foutname = tmpfile+"/tmpmerge.gz"
			fout = gzip.open(foutname, 'wt')
			
			merge_files(lst, fout)
			lsts.append(foutname)
			tmpfiles.append(foutname)
			fout.close()
	
	shutil.move(lsts[0], fnameout+"_perind.counts.gz")

def merge_files(fnames, fout):

	fopen = []
	for fname in fnames:
		if fname[-3:] == ".gz":
			fopen.append(gzip.open(fname))
		else:
			fopen.append(open(fname))

	finished = False
	N = 0
	while not finished:
		N += 1
		if N % 50000 == 0: 
			print(".")
		buf = []
		for f in fopen:
			ln = f.readline().split()
			if len(ln) == 0: 
				finished = True
				break
			chrom = ln[0]
			data = ln[1:]
			if len(buf) == 0:
				buf.append(chrom)
			buf += data
		if len(buf) > 0:
			if buf[0] == "chrom": 
				print("merging %d files"%(len(buf)-1))
			outStr = " ".join(map(lambda x: x.decode(), buf))
			fout.write(outStr + "\n")
		else:
			break

	print(" done.\n")
	for fin in fopen:
		fin.close()

#done
def cluster_intervals(E):
	''' Clusters intervals together. '''
 # E is a list of tuples<int, int>. the tuples are the coordinates of junctions. All of these junctions possess at least 3 reads, but are given without read counts here.
	E.sort()            # sorts clusters by start coord to optimise overlap - if you do not overlap a cluster starting at x, none at x+1 will overlap either
	if len(E) == 0: return [], []         # nothing to sort, return emprt
	current = E[0]
	Eclusters, cluster = [], []
	i = 0
	while i < len(E):            # iterate through list of clusters once
		if overlaps(E[i], current):       # check overlap between intron and 'current' item (start of intron and max of current end / all overlapping introns ends)
			cluster.append(E[i])       # overlapping intron is appended to the cluster
		else:              # since clusters are sorted by start coordinates, once one intron does not overlap, no other will, and we need to start a new cluster
			Eclusters.append(cluster)        # output current cluster (list of overlapping introns) to Eclusters object
			cluster = [E[i]]        # set cluster to start anew from current intron
		current = (E[i][0], max([current[1], E[i][1]]))  # current object is defined as the maximum length spanned by the start of the current item and the largest end coordinate of all overlapping introns
		i += 1

	if len(cluster) > 0: Eclusters.append(cluster)    # while loop will complete with last cluster unprocessed ('while' will stop at least intron) - if there is a current cluster it must be saved to the objest
	return Eclusters, E          # two possibly uneven items - the collated clusters and the original clusters, both lists. E is never used.

#done
def overlaps(A,B):
	'''
	Checks if A and B overlaps
	'''
	if A[1] < B[0] or B[1] < A[0]: return False # If end of first smaller than start of second, or end of second smaller than end of first, it's impossible for them to overlap. Else they do.
	else: return True         

#done
def refine_linked(cluster):
	''' Refines linked cluster (overlapping) by checking which of its introns share exact splice sites. 
	Input: list of introns. Output: list of cluster (which are lists of introns) that share splice sites.
	This may fragment existing cluster, which thus require validation, again '''
	#since we sorted before creating the file, these entries come in order of start position
	unassigned = [x for x in cluster[1:]] #ensures this is always a list
	current = [cluster[0]]
	splicesites = set([current[0][0][0],current[0][0][1]]) #start and end of the 'current' cluster first intron. Set, so only unique values remain.
	newClusters = []
	while len(unassigned) > 0:
		finished = False
		while not finished:
			finished = True
			torm = []
			for intron in unassigned:
				inter, count = intron
				start, end = inter
				if start in splicesites or end in splicesites: #rather than an overlap, it is now tested whether the intron shares splice sites with the... first intron plus any up to this point? I think??
					current.append(intron)
					splicesites.add(start)
					splicesites.add(end)
					finished = False
					torm.append(intron)
			for intron in torm:
				unassigned.remove(intron) #those introns matched to a refined cluster are taken from the list. this may skip several introns
		newClusters.append(current)
		current = []
		if len(unassigned) > 0:
			current = [unassigned[0]]		#grabs first unassigned intron, which may be one skipped when an earlier cluster was refined
			splicesites = set([current[0][0][0],current[0][0][1]]) #resets splice sites to match start and end of the current intron
			unassigned = unassigned[1:]
	return newClusters

#done
def refine_cluster(clu, cutoff, readcutoff, chr="NA"):
	''' for each exon in the cluster compute the ratio of reads, if smaller than cutoff, remove and recluster '''
	intervals = []
	reCLU = False
	totN = 0
	nIntronsFailProps = 0
	nIntronsFailCount = 0
	IsBIN1 = False
	if (chr == BIN1[0] and min(cl, key = lambda x: x[0])[0] >= BIN1[1] and max(cl, key=lambda x: x[1])[1] <= BIN1[2]): IsBIN1 = True
	for inter, count in clu: totN += count
	for inter, count in clu:
		testCounts = count >= readcutoff
		testProps  = count/float(totN) >= cutoff
		if (testCounts and testProps): #only admits refined clusters whose percentage of total counts is above cutoff and number of counts is above cutoff
			intervals.append(inter)
			dic[inter] = count
		else:
			if (not testCounts): nIntronsFailCount+=1
			if (not testProps): nIntronsFailProps +=1
			reCLU = True										# at least one intron was below one of the cutoff; which may collapse the entire assignment and requires a re-do.
	if ((not testCounts) or (not testProps)):
		logText = f"refine_cluster, cluster {chr}:{clu[0]} ({len(clu)} elements): "
		if (not testCounts): logText += f"{nIntronsFailCount} failed minimum counts. "
		if (not testProps):  logText += f"{nIntronsFailProps} failed to contain more than {cutoff*100}% of summed cluster reads. "
		logText += f"{len(intervals)} remaining introns will be reclustered."
		if (IsBIN1): logWrite(BIN1log, logText)
		logWrite(logFile, logText)
	
	#intervals now contains all introns that have passed the min count and min proportion tests
	if len(intervals) == 0: return []
	
	# This makes sure that after trimming, the clusters are still good
	Atmp, B = cluster_intervals(intervals)			#we re-cluster (overlap) the remaining introns.
	A = []
	for cl in Atmp:
		for c in refine_linked([(x,0) for x in cl]): # we refine linkage again, i.e. clusters with only exact overlap
			if len(c) > 0:							 # this does NOT filter single-intron clusters
				A.append([x[0] for x in c])
	
	if len(A) == 1:									# if this leaves only one cluster with 1+ introns:
		rc = [(x, dic[x]) for x in A[0]]			#takes all introns for this cluster? I think?
		if len(rc) > 1:								# if there's more than one intron, which... should be certain? Ah no it's 
			if reCLU:								# if any intron was removed because it violated min reads or min read percentage
				return refine_cluster([(x, dic[x]) for x in A[0]], cutoff, readcutoff, chr) # oh god recursion
			else:
				return [[(x, dic[x]) for x in A[0]]]
		else:
			return []
	NCs = []
	for c in A:
		if len(c) > 1:
			NC = refine_cluster([(x, dic[x]) for x in c], cutoff, readcutoff, chr)
			NCs += NC
	return NCs

if __name__ == "__main__":
	#options
	os.chdir(f"{os.path.dirname(os.path.abspath(__file__))}\\juncs")
	chromLst = ["chr%d"%x for x in range(1,23)]+['chrX','chrY']+["%d"%x for x in range(1,23)]+['X','Y'] 
	outPrefix = "HACKv1"
	rundir = "."
	maxIntronLen = 500000
	checkchrom = False
	useStrand = False
	minratio = 0.001
	minreads = 10
	runName = "%s/%s"%(rundir, outPrefix)

	def ts(): return datetime.datetime.now().strftime("[%y-%m-%d %H:%M:%S] --")
	def logWrite(logFile, logStr): 
		outStr = f"{ts()} {logStr}\n"
		sys.stdout.write(outStr)
		logFile.write(outStr)
	logFile = open(f"{rundir}/{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}_{outPrefix}.log", "w")
	BIN1log = open(f"{rundir}/{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}_{outPrefix}_BIN1.log", "w")
	
	logWrite(logFile, "Starting custom clustering run...")
	#juncfiles = ["Sample-161206-01_1W_Aligned.out.bam.junc", "Sample-161206-02_1W_Aligned.out.bam.junc", "Sample-161207-01_1W_Aligned.out.bam.junc", "Sample-161207-02_1W_Aligned.out.bam.junc", 
 #    "Sample-161212-01_1W_Aligned.out.bam.junc", "Sample-161208-01_1W_Aligned.out.bam.junc", "Sample-161208-02_1W_Aligned.out.bam.junc", "Sample-161206-08_AD_Aligned.out.bam.junc", 
 #    "Sample-161208-08_AD_Aligned.out.bam.junc", "Sample-161207-06_AD_Aligned.out.bam.junc", "Sample-161207-07_AD_Aligned.out.bam.junc", "Sample-170110-05_AD_Aligned.out.bam.junc", 
 #    "Sample-161207-10_AD_Aligned.out.bam.junc", "Sample-170124-07_AD_Aligned.out.bam.junc", "Sample-170221-01_AD_Aligned.out.bam.junc"] 
	juncfiles = ["Sample-161206-01_1W_Aligned.out.bam.junc", "Sample-161206-02_1W_Aligned.out.bam.junc"]
	libl = []
	for junc in juncfiles:
		junc = junc.strip()
		try: open(junc)
		except: 
			print("%s does not exist... check your junction files.\n"%junc)
			exit(0)
		libl.append(junc)
	# END OF ORIGINAL __main__
	logWrite(logFile, f"Found {len(libl)} .junc files")
	## start of pool_junc_reads()
	outFile = "%s/%s_pooled"%(rundir,outPrefix)
	by_chrom = {}
	BIN1Introns = 0
	for x in libl:
		lib = x.strip()
		if not os.path.isfile(lib): continue
		print("scanning %s...\n"%lib)

		if lib[-3:] == ".gz": F = gzip.open(lib)
		else: F = open(lib, "r")
		for ln in F: 
			lnsplit=ln.split()
			if len(lnsplit)<6: 
				print("Error in %s \n" % lib)
				continue
			chrom, A, B, dot, counts, strand = lnsplit					# line is split into chromosome, start, stop, a dot (because file format?) the counts supporting this intron and the strand
			if (chrom == BIN1[0] and int(A) <= BIN1[2] and int(B) >= BIN1[1] ): 
				BIN1Introns += 1
			if not useStrand: strand = "NA"
			if checkchrom and (chrom not in chromLst): continue         # removes all chrs that aren't (chr)+[1-23XY]
			A, B = int(A), int(B)+1										# A and B are integers now, B gets one added
			if B-A > int(maxIntronLen): continue						# introns over max length abandoned here
	     #Saving introns below max length to master list by_chrom
			try: by_chrom[(chrom,strand)][(A,B)] = int(counts) + by_chrom[(chrom, strand)][(A,B)]   # entry already exists, adding counts to pre-existing entry for interval A B
			except: 
				try: by_chrom[(chrom,strand)][(A,B)] = int(counts)      # IndexError; creating entry
				except: by_chrom[(chrom, strand)] = {(A,B):int(counts)} # KeyError: creates a tuple:int dictionary of (start,stop): supporting counts
	logWrite(BIN1log, f"Found {BIN1Introns} introns inside BIN1 across {len(libl)} files.")
	
	Ncluster = 0
	NCLTTT = 0
	BIN1Introns = []
	BIN1Clusters = []

	fout = open(outFile, 'w')

	Refined_Clusters = {}

	for chrom in by_chrom.keys():                    # this is per chromosome
		read_ks = [k for k,v in by_chrom[chrom].items() if v >= 3]          # takes all junctions on this chromosome that have at least 3 reads, whilst discarding counts, I think?
		read_ks.sort()
		

		if (chrom[0] == BIN1[0]):																	# blatantly assuming unstranded data here, sorry
			BIN1Introns = [ k for k in read_ks if k[0] <= BIN1[2] and k[1] >= BIN1[1] ]
			logWrite(BIN1log, f"{len(BIN1Introns)} BIN1 introns have more than 3 reads.")
		#else: continue																				# DEBUG HACK - bypasses rest of loop if chr != 2

		NIntrons = len(by_chrom[chrom])
		NCLTTT += ( NIntrons - len(read_ks) )
		logWrite(logFile, f"Chromosome {chrom}: Total {NIntrons} introns. Discarded {NIntrons - len(read_ks)} (less than 3 reads supporting intron). {len(read_ks)} introns going to overlap check.")
		
		clus = cluster_intervals(read_ks)[0]              # Sorta and collates clusters for this chromosome; returns ([collated clusters], [original clusters])
		
		if (chrom[0] == BIN1[0]):
			for clu in clus:
				if (min(clu, key = lambda x: x[0])[0] >= BIN1[1] and max(clu, key=lambda x: x[1])[1] <= BIN1[2]):
					BIN1Clusters.append(clu)
			logWrite(BIN1log, f"{len(BIN1Introns)} BIN1 introns have been sorted into {len(BIN1Clusters)} overlap cluster(s).")
		nClusters = 0
		for clu in clus: nClusters += len(clu)
		logWrite(logFile, f"\t{len(read_ks)} introns have been collated into {len(clu)} clusters. Sanity check: {nClusters} total introns within clusters.")
		nOneIntronClusters = 0
		nMultiIntronClusters = 0
		BIN1Clusters = []

		if not chrom in Refined_Clusters: Refined_Clusters[chrom]=[]
		
		Refined_Clusters[chrom].append(clu)
		for clu in clus:
			if len(clu) > 1: # if cluster has more than one intron - can't show alternative splicing if there's no alternative at these coords?
				if not chrom in Refined_Clusters: Refined_Clusters[chrom]=[]
				Refined_Clusters[chrom].append(clu)
				nMultiIntronClusters += 1
				buf = '%s:%s '%chrom                # reinserts chromosome information; collated cl 
				for interval, count in [(intron, by_chrom[chrom][intron]) for intron in clu]:	#for each cluster with multiple members... this grabs the cluster's data?
					buf += "%d:%d" % interval + ":%d"%count+ " "
					if (chrom[0] == BIN1[0] and interval[1] >= BIN1[1] and interval[0] <= BIN1[2]): 
						BIN1Clusters.append((interval, count))
				fout.write(buf+'\n')
			Ncluster += 1
		if (len(BIN1Clusters) > 0):
			logWrite(BIN1log, f"{len(BIN1Clusters)} BIN1 multi-intron overlap cluster(s) have been written to '_pooled' file")
			BIN1Clusters = []
		nOneIntronClusters = len(clu) - nMultiIntronClusters
		logWrite(logFile, f"\tDiscarded {nOneIntronClusters} single-intron clusters, leaving {nMultiIntronClusters} for differential usage analysis.")
	logWrite(logFile, f"Wrote {Ncluster} clusters total to {outFile}")
	logWrite(logFile, f"A total of {NCLTTT} introns were discarded (less than 3 reads supporting their existence).")
	BIN1log.flush()
	fout.flush()
	fout.close()
	#end of pool_junc_reads()

	del by_chrom
	del NCLTTT
	del Ncluster

	#start of refine_clusters()
	inFile = "%s/%s_pooled"%(rundir,outPrefix)
	outFile = "%s/%s_refined"%(rundir,outPrefix)

	fout = open(outFile,'w')
	Ncl = 0
	for ln in open(inFile):
		clu = []
		totN = 0
		chrom = ln.split()[0]
		for ex in ln.split()[1:]:
			A, B, N = ex.split(":")
			clus.append(((int(A),int(B)), int(N)))
			totN += int(N)
		if totN < minreads: continue
		
		clus = refine_linked(clu) #takes a single cluster within one chromosome; refine_linked returns list
		for clu in clus: #TODO check expected input, use Refined_Cluster design for below
			rclu = refine_cluster(clu, minratio, minreads)
			if len(rclu) > 0:
				for clu in rclu:	#each clu is a list of ((start, end), counts)
					buf = '%s ' % chrom
					for interval, count in clu:
						buf += "%d:%d" % interval + ":%d"%(count)+ " "
					Ncl += 1
					fout.write(buf+'\n')
	sys.stderr.write("Split into %s clusters...\n"%Ncl)
	fout.close()
	# end of refine_clusters() #

	Double_Refined_Clusters = []
	for chrom in Refined_Clusters.keys():
		chrom_Refined_Clusters = Refined_Clusters[chrom]
		for cluster in chrom_Refined_Clusters: #compare with clu
			for cl in refine_linked(cluster): #we call refind_linked() on all introns of current cluster. It checks for exact splice site overlaps and reclusters accordingly
				rc = refine_cluster(cl, minratio, minreads, chrom[0]) 
				#these are new clusters whose introns overlap at splice sites. We check they meet out minimum read requirement 
				#and that each intron has >0.1% of the cluster's reads
				if len(rc) > 0:
					Double_Refined_Clusters.append(rc)
					for clu2 in rc:
						for interval, count in clu2:
							if (chrom[0]==BIN1[0] and interval[0] <= BIN1[2] and interval[1] >= BIN1[1] ): 
								logWrite(BIN1log, "A BIN1 cluster has survived refinement. Saving to buffer.")
	
	# start of sort_junctions() #
	exons, cluExons = {}, {}
	cluN = 0
	for ln in open("%s/%s_refined"%(rundir,outPrefix), 'r'): #opening refined clusters back up, splitting AGAIN.... jesus, really?
		chrom = ln.split()[0]
		chrom = tuple(chrom.split(":")) #chromosome strand
		cluN += 1
		for exon in ln.split()[1:]:
			A, B, count = exon.split(":") #are these really exons
			if chrom not in exons: exons[chrom] = {}
			exons[chrom][(int(A),int(B))] = cluN
			if cluN not in cluExons: cluExons[cluN] = []
			cluExons[cluN].append((chrom, A, B))

	merges = {}
	#for ll in libl:
	#	lib=ll.rstrip()
	#	if not os.path.isfile(lib): continue
	#	libN = lib
	#	if libN not in merges: merges[libN] = []
	#	merges[libN].append(lib) #what are you -doing-

	for ll in map(lambda x: x.rstrip(), libl):
		if not os.path.isfile(ll): continue
		if ll not in merges: merges[ll] = []
		merges[ll].append(ll)
	
	fout_runlibs = open(runName+"_sortedlibs",'w')

	#replace merges with Double_Refined_Clusters, I guess? compare.
	for libN in merges:
		libName = "%s/%s"%(rundir,libN.split('/')[-1])
		by_chrom = {}
		foutName = libName+'.%s.sorted.gz'%(runName.split("/")[-1])

		fout_runlibs.write(foutName+'\n')

		print("Sorting %s..\n"%libN)
		if len(merges[libN]) > 1: print("merging %s...\n"%(" ".join(merges[libN])))
		else:
			pass
		
		fout = gzip.open(foutName,'wt')
		fout.write("chrom %s\n"%libN.split("/")[-1].split(".junc")[0]) #column header: chrom space source BAM file (i.e. sample)

		for lib in merges[libN]:
			if lib[-3:] == ".gz": F = gzip.open(lib)
			else: F = open(lib)
			for ln in F:
				lnsplit=ln.split()
				if len(lnsplit)<6: 
					print("Error in %s \n" % lib)
					continue
				chrom, start, end, dot, count, strand = ln.split()
				if not useStrand: strand = "NA"
				if checkchrom and (chrom not in chromLst): continue
				chrom = (chrom, strand)
				if chrom not in by_chrom: by_chrom[chrom] = {}							#ensures that there's a dict to populate
				intron = (int(start), int(end)+1)
				if intron in by_chrom[chrom]: by_chrom[chrom][intron] += int(count)
				else: by_chrom[chrom][intron] = int(count)
		for clu in cluExons:
			buf = []
			ks = cluExons[clu]
			ks.sort()
			tot = 0
			for exon in ks:
				chrom, start, end = exon
				start, end = int(start), int(end)
				if chrom not in by_chrom: pass
				elif (start,end) in by_chrom[chrom]: tot += by_chrom[chrom][(start,end)]
			for exon in ks:
				chrom, start, end = exon
				start, end = int(start), int(end)
				chromID, strand = chrom
				if chrom not in by_chrom: buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID, start, end, clu, strand, tot))
				elif (start,end) in by_chrom[chrom]: buf.append("%s:%d:%d:clu_%d_%s %d/%d\n"%(chromID, start, end, clu, strand, by_chrom[chrom][(start,end)], tot))
				else: buf.append("%s:%d:%d:clu_%d_%s 0/%d\n"%(chromID,start, end,clu,strand, tot))
			fout.write("".join(buf))
		fout.close()
	fout_runlibs.close()
	# end of sort_junctions() #

	# start of merge_junctions() #
	fnameout = "%s/%s"%(rundir,outPrefix)
	flist = "%s/%s_sortedlibs"%(rundir, outPrefix)
	
	lsts = []
	for ln in open(flist):
		lsts.append(ln.strip())
	logWrite(logFile, f"Merging {len(lsts)} junction files...")
	
	# Change 300 if max open file is < 300
	N = min([300, max([100, int(len(lsts)**(0.5))])]) #** 0.5 - a square root??
	N = min([300, max([ 100, int(math.sqrt(len(lsts))) ]) ])

	tmpfiles = []
	while len(lsts) > 1: 
		clst = []
		
		for i in range(0, math.ceil( (len(lsts)/N)+1 )): 
			lst = lsts[N*i:N*(i+1)]
			if len(lst) > 0:
				clst.append(lst)
		lsts = []													#splits full data into N chunks?
	
		for lst in clst:
			if len(lst) == 0: continue
			tmpfile = tempfile.mktemp()
			os.mkdir(tmpfile)
			foutname = tmpfile+"/tmpmerge.gz"
			fout = gzip.open(foutname,'wt')
			
			merge_files(lst, fout)
			lsts.append(foutname)
			tmpfiles.append(foutname)
			fout.close()
	
	shutil.move(lsts[0], fnameout+"_perind.counts.gz") #rename
	# end of merge_junctions() #

	# start of get_numers() # 
	fname = "%s/%s_perind.counts.gz"%(rundir,outPrefix)
	fnameout = "%s/%s_perind_numers.counts.gz"%(rundir,outPrefix)
	input_file=gzip.open(fname, 'r')
	fout = gzip.open(fnameout, 'wb')
	first_line=True

	for l in input_file:
		if first_line:
			fout.write(b" ".join(l.strip().split(b" ")[1:])+b"\n") # print the sample names
			first_line=False
		else:
			l=l.strip()
			words=l.split(b" ")
			fout.write(words[0]+ b" "+ b" ".join( [ g.split(b"/")[0] for g in words[1:] ] ) +b'\n')
	input_file.close()
	fout.close()
	# end of get_numers() #
	# end of main() #
	logFile.close()
	BIN1log.close()