using System;
using System.ComponentModel;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Forms;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using WinForms = System.Windows.Forms;

namespace FastQCRipper
{
	/// <summary>
	/// Interaction logic for MainWindow.xaml
	/// </summary>
	public partial class MainWindow : Window
	{
		const string NODIRSTR = "(No directory selected)";

		FolderBrowserDialog OpenDialog = new FolderBrowserDialog();
		SolidColorBrush SolidBrush = new SolidColorBrush();
		SolidColorBrush TransparentBrush = new SolidColorBrush();
		//BackgroundWorker ZIPRipper = new BackgroundWorker(); //TODO: ~multithreading~ (which probably isn't worth the effort here)

		public MainWindow(){
			InitializeComponent();
			SolidBrush.Opacity = 1.0;
			SolidBrush.Color = Color.FromRgb(0, 0, 0);
			TransparentBrush.Opacity = 0.5;
			TransparentBrush.Color = Color.FromRgb(0, 0, 0);
			TextBoxDestinationDirectory.Foreground = TransparentBrush;
			TextBoxSourceDirectory.Foreground = TransparentBrush;
			#if DEBUG
				TextBoxSourceDirectory.Text = @"D:\Users\User\Desktop\Cache";
				TextBoxDestinationDirectory.Text = @"D:\Users\User\Desktop\Cache";
			#endif
		}

		private void ButtomCompute_Click(object sender, RoutedEventArgs e){
			ButtomCompute.IsEnabled = false;
			ButtonSelectDestinationDirectory.IsEnabled = false;
			ButtonSelectSourceDirectory.IsEnabled = false;

			if ( TextBoxSourceDirectory.Text == NODIRSTR || TextBoxDestinationDirectory.Text == NODIRSTR ) {
				System.Windows.MessageBox.Show("Please select both a source and destination directory.");
				return;
			}
			
			DirectoryInfo SourceDirectory = new DirectoryInfo(TextBoxSourceDirectory.Text);
			if (!SourceDirectory.Exists) { System.Windows.MessageBox.Show(string.Format("Source directory \"{0}\" does not exist. Please use the \"Select...\" button to navigate to a real directory.", SourceDirectory.FullName)); return; }
			DirectoryInfo DestinationDirectory = new DirectoryInfo(TextBoxDestinationDirectory.Text);
			if (!DestinationDirectory.Exists) { System.Windows.MessageBox.Show(string.Format("Destination directory \"{0}\" does not exist. Please use the \"Select...\" button to navigate to a real directory.", SourceDirectory.FullName)); return; }

			ButtomCompute.Visibility = Visibility.Hidden;
			ButtonSelectDestinationDirectory.Visibility = Visibility.Hidden;
			ButtonSelectSourceDirectory.Visibility = Visibility.Hidden;
			TextBoxDestinationDirectory.Visibility = Visibility.Hidden;
			TextBoxSourceDirectory.Visibility = Visibility.Hidden;
			ProgressBarZIPRip.Visibility = Visibility.Visible;
			LabelProgress.Visibility = Visibility.Visible;

			FileInfo[] ZIPFiles = SourceDirectory.GetFiles("*_fastqc.zip");
			string[] FQZDataSplit = null;
			string rootname;
			List<string> Summary = new List<string>();
			int nFiles = ZIPFiles.Length;
			for (int i = 0; i < nFiles; i++) {
				LabelProgress.Content = string.Format("Processing file {0} of {1} ({2:0.0}%)...", i, nFiles, ((double)i/nFiles)*100);
				FileInfo FastQCZipFile = ZIPFiles[i];
				Summary.Add("Module Name\tStatus\n");
				//Extract the stats from their containing file
				using (ZipArchive FQZip = ZipFile.OpenRead(FastQCZipFile.FullName)) {
					rootname = FastQCZipFile.Name.Substring(0, FastQCZipFile.Name.Length - 4);
					ZipArchiveEntry FQZData = FQZip.Entries.Where(n => n.Name == "fastqc_data.txt").First(); //Filenames are without directory
					FQZDataSplit = new StreamReader(FQZData.Open()).ReadToEnd().Split(new string[] { ">>" }, StringSplitOptions.RemoveEmptyEntries);
				}
				//Probably a bad idea but hey it's only ~270KB so lets do this in-memory
				//File contains _all_ records, encapsulated by >>MODULE NAME\t<pass|fail>\n[...]>>END MODULE, thus, split by '>>'
				FQZDataSplit = FQZDataSplit.Where(n => (n != "END_MODULE\n")).Skip(1).ToArray();
				//now each line represents one dataset, to be written to one file
				//(first entry is the fast qc version info, and thus skipped)
				foreach (string FQZDataFile in FQZDataSplit) {
					string[] FQZSplit = FQZDataFile.Split('\n');
					string moduleName = FQZSplit[0]; //first line contains the module name and pass/fail
					Summary.Add(moduleName);
					moduleName = rootname + "_" + moduleName.Split('\t')[0] + ".txt"; //Create output file named based on the input file and current module
					FQZSplit[1] = FQZSplit[1].Replace("#", string.Empty); //correct headers
					try {
						FileInfo OutputFile = new FileInfo(System.IO.Path.Combine(DestinationDirectory.FullName, moduleName));
						using (StreamWriter OutputFileWriter = new StreamWriter(OutputFile.Create())) { OutputFileWriter.Write(string.Join("\n", FQZSplit.Skip(1))); } //dump stats to file
					}
					catch (UnauthorizedAccessException) { System.Windows.MessageBox.Show("The program is unable to gain access to the target directory. Please select a new output directory or re-run the application with elevated privileges."); return; }
					catch (PathTooLongException) { System.Windows.MessageBox.Show("The program cannot create the output file, as its full path would be longer than the system maximum. Please select a different output folder with a shorter path."); return; }
				}

				//repeat data dump for summary
				try {
					FileInfo OutputFile = new FileInfo(System.IO.Path.Combine(DestinationDirectory.FullName, (rootname + "_Summary.txt")));
					using (StreamWriter OutputFileWriter = new StreamWriter(OutputFile.Create())) { foreach (string s in Summary) { OutputFileWriter.Write(s + "\n"); } }
				}
				catch (UnauthorizedAccessException) { System.Windows.MessageBox.Show("The program is unable to gain access to the target directory. Please select a new output directory or re-run the application with elevated privileges."); }
				catch (PathTooLongException) { System.Windows.MessageBox.Show("The program cannot create the output file, as its full path would be longer than the system maximum. Please select a different output folder with a shorter path."); }
				//Update GUI - supposed to be a BackgroundWorkerThread...
			}

			ButtomCompute.Visibility = Visibility.Visible;
			ButtonSelectDestinationDirectory.Visibility = Visibility.Visible;
			ButtonSelectSourceDirectory.Visibility = Visibility.Visible;
			TextBoxDestinationDirectory.Visibility = Visibility.Visible;
			TextBoxSourceDirectory.Visibility = Visibility.Visible;
			ProgressBarZIPRip.Visibility = Visibility.Hidden;
			LabelProgress.Visibility = Visibility.Hidden;
			ButtomCompute.IsEnabled = true;
			ButtonSelectDestinationDirectory.IsEnabled = true;
			ButtonSelectSourceDirectory.IsEnabled = true;
			System.Windows.MessageBox.Show("Operation complete. Extracted " + nFiles + " FastQC files into individual .txt files.\nHave a nice day!");
			//Done!
		}

		private void ButtonSelectSourceDirectory_Click(object sender, RoutedEventArgs e){ DoFolderBrowse(TextBoxSourceDirectory); }

		private void ButtonSelectDestinationDirectory_Click(object sender, RoutedEventArgs e){ DoFolderBrowse(TextBoxDestinationDirectory); }

		private void DoFolderBrowse(System.Windows.Controls.TextBox Destination) {
			WinForms.DialogResult dr = OpenDialog.ShowDialog();
			switch (dr){
				case WinForms.DialogResult.OK:
				case WinForms.DialogResult.Yes:
					Destination.Text = OpenDialog.SelectedPath;
					Destination.Foreground = SolidBrush;
					break;
				default:
					Destination.Text = NODIRSTR;
					Destination.Foreground = TransparentBrush;
					break;
			}
		}
	}
}

