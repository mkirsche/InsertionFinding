/*
 * Takes in a list of how many times each kmer in a set of insertions
 * occurs in both the reference and a set of reads, and outputs for
 * different thresholds t how many insertions have at least one kmer which
 * occurs 0 times in the reference and >= t times in the readset.
 *
 * The output is one line per threshold of the form: threshold num_detected
 */
import java.util.*;
import java.io.*;
public class ParseOutput {
public static void main(String[] args) throws IOException
{
	String fn = args[0];
	PrintWriter out = new PrintWriter(new File(fn + ".signal"));
	for(int t = 1; t<=80; t++)
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
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
		System.out.println("Threshold for signal = " + t);
		System.out.println("Number of SVs = " + tot);
		System.out.println("Number with signal = " + count);
		System.out.println();
		out.println(t+" "+count);
	}
	out.close();
}
}
