using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Code generic key/value parser for the KEY:VALUE pairs (See the GTFParser you did)

namespace FASTAToDataTable {
	class Program {
		static void Main(string[] args) {
#if DEBUG
			args = new string[] { @"D:\Users\User\Documents\Actual Documents\University of Manchester\PhD\Ovis_aries.Oar_v3.1.cdna.all.fa" };
#else
			if (args.Length != 1) { printUsage(); return; }
#endif
			FileInfo FASTAFile = new FileInfo(args[0]);
			if (!FASTAFile.Exists) { Console.Out.Write("Error: Cannot open '" + args[0] + "': File does not exist."); return; }
			string Outfile = FASTAFile.DirectoryName + @"\" + DateTime.Now.ToString("yyMMdd-HHmm") + "_RFormat-FASTA-" + FASTAFile.Name + ".txt";
			try {
				using (Stream GTFFileStream = File.OpenRead(FASTAFile.FullName))
				using (Stream OutFileStream = File.OpenWrite(Outfile)) {
					StreamReader FastaReader = new StreamReader(GTFFileStream);
					StreamWriter OutfileWriter = new StreamWriter(OutFileStream);
					string currentLine;
					string currentSequence = "";
					OutfileWriter.Write("ParentID\tTranscriptID\tStarts/Stops\tStrand\tSequence" + Environment.NewLine);
					while ((currentLine = FastaReader.ReadLine()) != null) {
						if (currentLine[0] == '>') { //>0 MSTRG.1.1 1. 1-795 <- example
							if (currentSequence != "") {
								OutfileWriter.Write(currentSequence + Environment.NewLine); //Data.Table in R will NOT read inconsistent line endings.
								OutfileWriter.Flush();
								currentSequence = "";
							}
							string[] HeaderData = currentLine.TrimStart(new char[] { '>', ' ', '\t' }).Split(' ');
							if (HeaderData[1] == "cdna") {
								string TranscriptID = HeaderData[0].Substring(0, 3) == "ENS" ? HeaderData[0] : string.Join(".", HeaderData[0].Split('.').Take(HeaderData[0].Split('.').Length - 1));
								string GeneID = HeaderData.Where(x => x.Split(':')[0] == "gene").ToArray()[0];
								//string Location = HeaderData.Where(x => x.Split(':')[0] == "chromosome").ToArray()[0]; //missing in some cases, I guess? or wrong field. 
								OutfileWriter.Write(string.Format("{0}\t{1}\t\t\t", TranscriptID, GeneID));
							}
							else {
								string ParentID = HeaderData[1].Substring(0, 3) == "ENS" ? HeaderData[1] : string.Join(".", HeaderData[1].Split('.').Take(HeaderData[1].Split('.').Length - 1));
								OutfileWriter.Write(string.Format("{0}\t{1}\t{2}\t{3}\t", ParentID, HeaderData[1], HeaderData[3], HeaderData[2].TrimStart('1')));
							}
						}
						else { currentSequence += currentLine.Trim(); }
					}//while currentLine != null
				}
			}//try
			catch (Exception e) {
				Console.WriteLine("Error: " + e.Message + "\n" + e.Source);
			}//catch
		} //static void Main

		private static void printUsage()
		{
			Console.Out.WriteLine("Converts a .fa file into a data-table friendly format. Usage: FASTAToDataTable.exe <file.fa>");
		}
	}// class Program
}
