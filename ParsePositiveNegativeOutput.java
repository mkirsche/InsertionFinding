/*
 * Takes in a list of how many times each kmer in two sets of insertions
 * (a positive set and a negative set) occurs in both the reference
 * and a set of reads, and outputs for different thresholds t the
 * precision and recall (assuming an example is called as positive if it 
 * has at least one kmer not occurring in the reference but occurring >= t
 * times in the set of reads).
 *
 * This takes in two files - one for the positive examples, and another for
 * the negative examples
 *
 * The output is one line per threshold of the form: precision recall
 */
import java.util.*;
import java.io.*;
public class ParsePositiveNegativeOutput {
public static void main(String[] args) throws IOException
{
	String fnPositive = args[0];
	String fnNegative = args[1];
	String outFile = args[2];
	String[] fns = new String[]{fnPositive, fnNegative};
	PrintWriter out = new PrintWriter(new File(outFile));
	for(int t = 1; t<=80; t++)
	{
	    int[] tots = new int[2];
	    int[] counts = new int[2];
	    for(int file = 0; file < 2; file++)
	    {
		    Scanner input = new Scanner(new FileInputStream(new File(fns[file])));
		    int tot = 0;
		    int count = 0;
		    while(input.hasNext())
		    {
			    String line = input.nextLine();
			    if(line.startsWith("Adding")) continue;
			    if(line.charAt(0) < 'A' || line.charAt(0) > 'T') continue;
			    tot++;
			    String seq = line;
			    String[] gFreq = input.nextLine().split(" ");
			    String[] rFreq = input.nextLine().split(" ");
			    input.nextLine();
			    int n = gFreq.length;
			    int[] gs = new int[n], rs = new int[n];
			    for(int i = 0; i<n; i++)
			    {
				    gs[i] = Integer.parseInt(gFreq[i]);
				    rs[i] = Integer.parseInt(rFreq[i]);
			    }
			    boolean signal = false;
			    for(int i = 0; i<n; i++)
			    {
				    if(gs[i] == 0 && rs[i] >= t) signal = true;
			    }
			    if(signal) count++;
		    }
		    tots[file] = tot;
		    counts[file] = count;
		}
		System.out.println("Threshold for signal = " + t);
		System.out.println("Number of SVs = " + (tots[0] + tots[1]));
		int tp = counts[0];
		int fn = tots[0] - counts[0];
		int tn = tots[1] - counts[1];
		int fp = tots[1];
		double precision = 1.0 * tp / (tp + fp);
		double recall = 1.0 * tp / (tp + fn);
		System.out.println("Precision = " + precision);
		System.out.println("Recall = " + recall);
		System.out.println();
		out.println(precision+" "+recall);
	}
	out.close();
}
}
