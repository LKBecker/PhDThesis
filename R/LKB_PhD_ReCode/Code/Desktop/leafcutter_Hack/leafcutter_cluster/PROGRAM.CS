using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace leafcutter_cluster {
	public static class Config {
		public const int READ_NOISE_THRESHOLD = 3;
		public const int MIN_CLUSTER_INTRON_COUNTS = 10;
		public const double MIN_INTRON_CLUSTER_PERCENT = 0.001;
		public const int MAX_INTRON_LENGTS = 100000;
	}

	[DebuggerDisplay("Chr{Chr}[{Strand}]:{Start}-{End} Cluster: {Count} introns represent {TotalCounts} total counts")]
	class Cluster : HashSet<Intron> {
		internal string Chr;
		internal string Strand;
		internal bool isReducedToSites = false;
		internal readonly bool trackSize = false;
		internal readonly bool trackChr = false;

		internal int TotalCounts { get { return this.Sum( x => x.TotalCounts); } }
		public int Start = int.MaxValue;
		public int End = 0;

		public new void Add(Intron item) {
			if (trackChr && item.Chr != this.Chr) {
				throw new ArgumentException(string.Format("Cannot add to cluster: cluster is single-chromosome {0} but intron is from chromosome {1}", this.Chr, item.Chr));
			}
			if (this.Contains(item)) { //note that Intron.IsEqual does NOT check counts! This is by design.
				TryGetValue(item, out Intron currentValue);
				//#if DEBUG
				//				if (item.Chr == "2" && item.Strand=="+" && item.Start >= 117267199 && item.End <= 117314475) {
				//					Console.WriteLine("Uniting two BIN1 introns at {0}. Dict counts: {1}, incoming counts {2}", item, currentValue.Counts, item.Counts);
				//				}
				//#endif
				foreach (KeyValuePair<string, int> kv in item.Counts) { currentValue.AddCounts(kv.Key, kv.Value); } 
				//we already have an intron at these coordinates; unifying counts
				//this.Remove(item);
				//base.Add(currentValue); //using the underlying HashSet<Intron>'s Add() function to avoid recursion.
			} else {
				base.Add(item);
				//instead of a LINQ query across all elements each time we check Start or End (e.g. on each overlap) we do it on each element added. 
				//TODO is this less efficient to construct, more efficient to run Overlap()?
				if (trackSize) {
					if (item.Start < this.Start) { this.Start = item.Start; }
					if (item.End > this.End) { this.End = item.End; }
				}
			}
		}

		//constructors
		public Cluster(string chr, string strand) {
			this.Chr = chr;
			this.Strand = strand;
		}

		public Cluster(string chr, string strand, bool trackSize) {
			this.Chr = chr;
			this.Strand = strand;
			this.trackSize = trackSize;
		}

		internal List<Cluster> ClusterBySpliceSites() {
			//carries out the same function as refind_linked(cluster), except for the recluster part needing a different algo
			HashSet<int> MultiIntronSites = this.Select(x => x.Start).Concat(this.Select(y => y.End)).GroupBy(z => z)
				.Select(a => new { Coord = a.Key, Count = a.Count() }).Where(x => x.Count > 1).Select(x => x.Coord).ToHashSet();
			
			//For all positions where more than 1 intron(s) start or end, get those positions
			//we intersect with a subset of ourselves containing only overlappers, leaving only the overlappers 
			//(HashSet doesn't have an AddAll method and we can't replace our insides from inside the class)
			this.IntersectWith(this.Where(x => MultiIntronSites.Contains(x.Start) || MultiIntronSites.Contains(x.End)));
			if (this.Count < 2) { return new List<Cluster>(); }
			//List<Cluster> Algo2 = Program.MakeSpliceClustersVLeafcutter(this);
			//List<Cluster> Algo3 = Program.MakeSpliceClustersV2(this);
			return Program.MakeSpliceClustersV2(this);

			//TODO is this bad from a GC perspective
		}

		internal bool ApplyAndCheckIntronMinPercentage() {
			//TODO should this fail loudly or quietly call ReduceToSpliceSites)
			int tmpTotalCounts = this.TotalCounts;
			int tmpCounts = this.Count;
			//Here, we apply the same logic as in 
			List<Intron> GoodIntrons = this.Where(x => ((double)x.TotalCounts / tmpTotalCounts) >= Config.MIN_INTRON_CLUSTER_PERCENT && x.TotalCounts >= Config.MIN_CLUSTER_INTRON_COUNTS).ToList();
			if (GoodIntrons.Count < this.Count) { Program.nIntronsClusterPercentTooLow++; }
			IntersectWith(GoodIntrons);
			if (tmpCounts == this.Count) { return true; }
			// Program.Log(string.Format("{0} of {1} introns fail the minimum count % test. Removing and reclustering remains, if any...", this.Count, tmpCounts));
			return false;
		}

		public override string ToString() { return $"Chromosome {Chr}:{Strand} Cluster, {Count} introns"; }

	}

	[DebuggerDisplay("Chr{Chr}:{Start}-{End}[{Strand}]:{TotalCounts}")]
	class Intron : IComparable<Intron>{
		public readonly string Chr;
		public readonly int Start;
		public readonly int End;
		public readonly string Strand;
		public Dictionary<string, int> Counts;
		public int TotalCounts = 0;

		public Intron(int start, int end, string chr, int counts, string File, string strand = "NA"){
			this.Start  = start;
			this.End	= end + 1; //y tho
			if (this.End <= this.Start) { throw new ArgumentException("An introns' end must come at least 1 base after its start"); }
			this.Chr	= chr;
			Counts = new Dictionary<string, int>();
			this.Counts[File] = counts;
			this.TotalCounts = counts;
			this.Strand = strand;
		}

		public void AddCounts(string File, int counts) {
			if (Counts.ContainsKey(File)) {
				throw new ArgumentException("Junc file " + File + " should not contain counts for the same intron twice. Check input file(s).");
			}
			Counts[File] = counts;
			TotalCounts += counts;
		}

		internal bool Overlaps(Intron i){
			if (i.Chr != this.Chr || i.Strand != this.Strand) { return false; }
			if (this.Start < i.End && this.End > i.Start) { return true; }
			return false;
		}

		internal bool Overlaps(Cluster c) {
			if (c.Strand == "NA") {
				return c.Chr == this.Chr && c.Start <= this.End && c.End >= this.Start;
			} else {
				throw new NotImplementedException();
			}
		}

		//All equality operators for Intron ignore Counts!!
		public override bool Equals(object obj) {
			if (!(obj is Intron i)) { return false; }
			return i.Strand == this.Strand && i.Chr == this.Chr && i.Start == this.Start && i.End == this.End;
		}

		public override int GetHashCode() { return this.Chr.GetHashCode() ^ this.Strand.GetHashCode() ^ this.Start.GetHashCode() ^ this.End.GetHashCode(); }
		public override string ToString() { return string.Format("{0}:{1}-{2}[{3}]", this.Chr, this.Start, this.End, this.Strand); }

		public int CompareTo(Intron i) {  return this.Start.CompareTo(i.Start); } //used for Sorting
	}

	static class Program {
		public static int nIntronsTooLong = 0;
		public static int nIntronsClusterPercentTooLow = 0;
		public static int nSingleIntronClusters = 0;
		public static int nIntronsCombinedReadsTooLow = 0;

		public static void Log(string s) => Console.WriteLine(string.Format("{0} -- {1}", DateTime.Now.ToString("u"), s));

		public static void Log(string s, StreamWriter sw) {
			s = string.Format("{0} -- {1}", DateTime.Now.ToString("u"), s);
			Console.WriteLine(s);
			sw.WriteLine(s);
		}

		public static List<Cluster> MakeClustersByOverlap(List<Intron> Introns) {
			//TODO make TEST3, TEST4 with a few dummies and numbers that end in 5 for junction and 10 20 30 40 50 for counts
			List<Cluster> Clusters = new List<Cluster>();
			if (Introns.Count < 2) { return Clusters; }                             //can't make at least one overlap with less than 2 introns

			HashSet<string> Chromosomes = Introns.Select(x => x.Chr).ToHashSet();
			if (Chromosomes.Count > 1) { throw new ArgumentException("Cannot overlap introns: introns from more than one chromosome present in SortedSet!"); }

			Cluster c = new Cluster(Introns.First().Chr, "NA", true) { Introns.First() };   //creates new Cluster c and populates with current intron
			Intron[] remaining = Introns.Skip(1).ToArray();                                 
			for (int i = 0; i < remaining.Length; i++) {
				if (remaining[i].Overlaps(c)) { c.Add(remaining[i]); }  //Adding an intron to the cluster changes the cluster's Start and Stop, 
																		//so the cluster expands and its overlap does as well, iteratively
				else {
					if (c.Count > 1) {                                                      //OOPS WE'RE NOT SUPPOSED TO APPLY MIN COUNT JUST YET
						Clusters.Add(c);                                                    //We append our cluster to the list!!
					}
					c = new Cluster(remaining[i].Chr, "NA", true) { remaining[i] }; //Make new cluster, loaded with the first element that didn't overlap
				} //assuming these are sorted by Start, once one intron's start > end of current cluster, all subsequent starts will be, too, so no need to look at any element twice
			}
			Clusters.Add(c); //the final iteration of the loop won't get here (no object not to overlap with), so the last cluster must be added manually.
			return Clusters;
		}

        //C# implementation of leafcutters' clustering algorithm. Recreating this code taught me a lot about loops and how to retain my sanity.
		public static List<Cluster> MakeSpliceClustersV2(Cluster Introns) {
			List<Intron> RemainingIntrons = Introns.Skip(1).ToList();
			Cluster Current = new Cluster(Introns.First().Chr, "NA", true) { Introns.First() };
			HashSet<int> SpliceSites = new HashSet<int>() { Current.Start, Current.End };
			List<Cluster> Output = new List<Cluster>();

			while (RemainingIntrons.Count > 0) {
				bool NothingAdded = false;
				while (!NothingAdded) {
					NothingAdded = true;
					if (RemainingIntrons.Any(x => SpliceSites.Contains(x.Start) || SpliceSites.Contains(x.End))) {
						NothingAdded = false;
						List<Intron> Matches = RemainingIntrons.Where(x => SpliceSites.Contains(x.Start) || SpliceSites.Contains(x.End)).ToList();
						Matches.ForEach(x => { Current.Add(x); SpliceSites.Add(x.Start); SpliceSites.Add(x.End); });
						RemainingIntrons = RemainingIntrons.Where(x => !Matches.Contains(x)).ToList();
					}
				}
				Output.Add(Current);
				if (RemainingIntrons.Count > 0) {
					Current = new Cluster(RemainingIntrons.First().Chr, "NA", true) { RemainingIntrons.First() };
					RemainingIntrons = RemainingIntrons.Skip(1).ToList();
					SpliceSites = new HashSet<int>() { Current.Start, Current.End };
				}
			}
			return Output;
		}

		public static List<Cluster> MakeSpliceClustersVLeafcutter(Cluster Introns) {
			//Currently returns 3 clusters for first on Chr1, not... 1; use customer version above please!
			List<Intron> RemainingIntrons = Introns.Skip(1).ToList();
			Cluster Current = new Cluster(Introns.First().Chr, "NA", true) { Introns.First() };
			HashSet<int> SpliceSites = new HashSet<int>() { Current.Start, Current.End };
			List<Cluster> Output = new List<Cluster>();

			while (RemainingIntrons.Count > 0) {
				bool NothingAdded = false;
				while (!NothingAdded) {
					NothingAdded = true;							//We assume we won't find anything
					List<Intron> ToRemove = new List<Intron>();
					foreach (Intron intron in RemainingIntrons) {
						if (SpliceSites.Contains(intron.Start) || SpliceSites.Contains(intron.End)) {	//Each intron is checked for overlap
							Current.Add(intron);
							SpliceSites.Add(intron.Start);
							SpliceSites.Add(intron.End);
							NothingAdded = false;														//Since we -did- find something, all other introns must be evaluated 
							ToRemove.Add(intron);														//Else a far-aprt splice sites can overshadow splices inside it
						}
					}
					RemainingIntrons = RemainingIntrons.Where(x => !ToRemove.Contains(x)).ToList();
				}
				//the while !NotihngAdded loop has concluded - all introns were cycled through until none of them match any more; any left in RemainingIntorns MUST be a new cluster. so:
				//Add current
				Output.Add(Current); 
				//And if there's anything left, start afresh!
				if (RemainingIntrons.Count > 0) {
					Current = new Cluster(RemainingIntrons.First().Chr, "NA", true) { RemainingIntrons.First() };
					SpliceSites = new HashSet<int>() { Current.Start, Current.End };
					RemainingIntrons = RemainingIntrons.Skip(1).ToList();
				}
			}
			return Output;
		}

		static int Main(string[] args){
#if DEBUG
			//string ExePath = @"E:\Users\User\Documents\Actual Documents\University of Manchester\PhD\Code\R\LKB_PhD_ReCode\Code\Desktop\leafcutter_Hack\juncs-BHB";
			string ExePath = @"C:\PhD\Code\R\LKB_PhD_ReCode\Code\Desktop\leafcutter_Hack\juncs";
			//string FileName = Path.Combine(ExePath, "clusters.bin");
#else
			string ExePath = System.IO.Directory.GetCurrentDirectory();
#endif
			FileInfo[] juncFiles = new DirectoryInfo(ExePath).GetFiles().Where(n => n.Extension == ".junc").ToArray(); //load files generated using bam2junc.sh from original Leafcutter package
			//juncFiles = juncFiles.Where(n => n.Name.Contains("AD") | n.Name.Contains("1W")).ToArray();
			Cluster AllIntrons = new Cluster("NA", "NA"); //Collector for all entries from all files for initial stage, doesn't track size for slightly better performance
			//Has to be a Cluster instance, to merge counts from the same intron from multiple samples

			//open all junc files, read in all introns from them
			foreach (FileInfo juncFile in juncFiles) {
				using (FileStream fs = File.Open(juncFile.FullName, FileMode.Open))
				using (StreamReader sr = new StreamReader(fs)){
					int nIntronsRead = 0;
					while (!sr.EndOfStream) {
						string line = sr.ReadLine(); //example: 1	168572269	168579986	.	52	+
						if (line != null) {
							string[] lineSpl = line.Split('\t');
							int intron_start = int.Parse(lineSpl[1]);
							int intron_end = int.Parse(lineSpl[2]);
							if ((intron_end - intron_start) > Config.MAX_INTRON_LENGTS) {
								nIntronsTooLong++;
								continue;
							}
							if (lineSpl.Length == 6) { //we have strand information
								AllIntrons.Add(new Intron(intron_start, intron_end, lineSpl[0], int.Parse(lineSpl[4]), Path.GetFileNameWithoutExtension(juncFile.Name), lineSpl[5]));
							} else {
								AllIntrons.Add(new Intron(intron_start, intron_end, lineSpl[0], int.Parse(lineSpl[4]), Path.GetFileNameWithoutExtension(juncFile.Name), strand: "NA")); 
							}
							++nIntronsRead;
						} //line != null
					} //!sr.EndOfStream
					Log(String.Format("{0} total introns loaded from '{1}'", nIntronsRead, juncFile.Name));
				} //StreamReader sr, FileStream fs
			} //foreach juncFile
			
			//filter out introns with less than 3 reads
			int OldCount = AllIntrons.Count;
			IEnumerable<Intron> IntFiltTest = AllIntrons.Where(x => (x.TotalCounts >= Config.READ_NOISE_THRESHOLD));
			nIntronsCombinedReadsTooLow = OldCount - IntFiltTest.Count();
			Log(string.Format("{0} of {1} unique introns pass minimum count threshold (counts > {2})", IntFiltTest.Count(), OldCount, Config.READ_NOISE_THRESHOLD));
			AllIntrons = null;


            //EARLY TABLE DUMP
            string[] PossibleFiles = juncFiles.Select(n => Path.GetFileNameWithoutExtension(n.Name)).ToArray();
			StreamWriter countsWriter = new StreamWriter(Path.Combine(ExePath + "_perIntronCounts-Raw.out"));
			countsWriter.WriteLine("Intron.Chr\tIntron.Start\tIntron.Stop\t" + string.Join('\t', PossibleFiles));
			foreach (Intron iFinal in IntFiltTest) {
				if (iFinal.TotalCounts <= Config.MIN_CLUSTER_INTRON_COUNTS) { continue; }
				countsWriter.Write($"{iFinal.Chr}\t{iFinal.Start}\t{iFinal.End}");
				foreach (string PossibleFile in PossibleFiles) {
					if (iFinal.Counts.TryGetValue(PossibleFile, out int i)) { countsWriter.Write("\t" + i); }
					else { countsWriter.Write("\t0"); }
				}
				countsWriter.Write("\n");
			}
			countsWriter.Close();
            Log("Exported to <SAMPLE>-perIntronCounts-Raw.out file. Program stage 1 complete.");

            Environment.Exit(0); //this skips clustering and all that leafcutter-y stuff

			//Separate by chromosome, then sort by coords
			Dictionary<string, List<Cluster>> perChromClusters = new Dictionary<string, List<Cluster>>(); //perChromClusters["2"].Where( x => x.Start <= 117314475 && x.End >= 117267199)
			string[] Chromosomes = IntFiltTest.Select(x => x.Chr).Distinct().ToArray();
			//TODO string[] strands = AllIntrons.Select(x => x.Chr).Distinct().ToArray(); 
			//TODO chrCheck
			foreach (string Chromosome in Chromosomes) {
				//grab all introns that have same chr; order by Start coord; we use a SortedSet to ensure this doesn't mess up.
				//TODO sort once, use list? Should preserve insertion order
				List<Intron> chrIntrons = IntFiltTest.Where(x => x.Chr == Chromosome).OrderBy(x => x.Start).ThenBy( x=> x.End).ToList();
				if (chrIntrons.Count > 0) {
					List<Cluster> chrRoughClusters = MakeClustersByOverlap(chrIntrons);
					chrRoughClusters = chrRoughClusters.Where(x => x.Count > 1).ToList();
					nSingleIntronClusters += chrRoughClusters.Sum(x => x.Count) - chrIntrons.Count;
					perChromClusters[Chromosome] = chrRoughClusters;
				}
			} //foreach Chromosome
			Log(String.Format("Finished overlap. Final Result: {0} clusters", perChromClusters.Values.Sum( x => x.Count)));

			// end of pool_junc_reads()
			// ---- EVERYTHING ABOVE THIS LINE GIVES RESULTS IDENTICAL TO THOSE OF leafcutter_cluster.py ----

			//TODO: separate for now, eventually merge
			//StreamWriter countsWriter = new StreamWriter(Path.Combine(ExePath + "_perIntronCounts.out"));
			countsWriter = new StreamWriter(Path.Combine(ExePath + "_perIntronCounts-Refined.out"));
			countsWriter.WriteLine("Intron " + string.Join(' ', PossibleFiles));
			foreach (string Chromosome in Chromosomes) {
				Queue<Cluster> ClustersToRefine = new Queue<Cluster>(perChromClusters[Chromosome]);
				List<Cluster> chrRefinedClusters = new List<Cluster>();
				while (ClustersToRefine.Count > 0) {
					if (ClustersToRefine.TryDequeue(out Cluster c)) {		//Grab overlap cluster
						foreach (Cluster c2 in c.ClusterBySpliceSites()) {	//Reduce to splice site clusters (reduce to splice, recluster)
							if (!c2.ApplyAndCheckIntronMinPercentage()) {	//If any splice cluster failt Intron number or Intron percentage:
								ClustersToRefine.Enqueue(c2);               //Reinsert into Queue, which will recluster, then recheck (low-read/% introns are gone)
							} else {
								//chrRefinedClusters.Add(c2);
								foreach (Intron iFinal in c2) {
									countsWriter.Write($"{iFinal.Chr}:{iFinal.Start}:{iFinal.End}");
									foreach (string PossibleFile in PossibleFiles) {
										if (iFinal.Counts.TryGetValue(PossibleFile, out int i)) { countsWriter.Write(" " + i); }
									}
									countsWriter.Write("\n");
								}
							}
						}
						countsWriter.Flush();
					} //else { break; }
				}
				chrRefinedClusters = chrRefinedClusters.OrderBy(x => x.Start).ThenBy(y => y.End).ToList();
				perChromClusters[Chromosome] = chrRefinedClusters;
			}
			countsWriter.Close();
			//Now we should have a Dictionary of < chromosome, List<Clusters> for splice junctions. 
			//We can now use this to extract the relevant counts from the original .junc files, and make a per-sample table. I guess.
			Dictionary<string, List<Cluster>> perFileSpliceCounts = new Dictionary<string, List<Cluster>>();
			//Final format: chrom <sample 1..n>
			//				1		
			return 0;
		}
	}
}
