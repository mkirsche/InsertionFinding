/*
 * Takes in the positive and negative examples, as well as a threshold for calling
 * Outputs the calls within each set of examples as well as additional information
 * about the sequences which make up insertion calls
 */
import java.util.*;
import java.io.*;
public class GetPositiveExamples {
public static void main(String[] args) throws IOException
{
	String fn1 = args[0];
	String fn2 = args[1];
	String[] fns = new String[]{fn1, fn2};
	ArrayList<String> seqs1 = getAllSeqs(fn1), seqs2 = getAllSeqs(fn2), allseqs = new ArrayList<String>();
	allseqs.addAll(seqs1);
	allseqs.addAll(seqs2);
	int t = Integer.parseInt(args[2]);
	HashMap<String, Integer> freqMap = null;
	PrintWriter out1 = new PrintWriter(new File(fn1 + ".calls." + t));
	PrintWriter out2 = new PrintWriter(new File(fn2 + ".calls." + t));
	Scanner input1 = new Scanner(new FileInputStream(new File(fn1)));
	Scanner input2 = new Scanner(new FileInputStream(new File(fn2)));
	Scanner[] inputs = new Scanner[]{input1, input2};
	PrintWriter[] outs = new PrintWriter[]{out1, out2};
	int offset = 0;
	for(int p = 0; p<inputs.length; p++)
	{
		int tot = 0;
	    int count = 0;
	    while(inputs[p].hasNext())
	    {
		    String line = inputs[p].nextLine();
		    if(line.startsWith("Adding")) continue;
		    if(line.charAt(0) < 'A' || line.charAt(0) > 'T') continue;
		    tot++;
		    String seq = line;
		    String[] gFreq = inputs[p].nextLine().split(" ");
		    String[] rFreq = inputs[p].nextLine().split(" ");
		    inputs[p].nextLine();
		    int n = gFreq.length;
		    int[] gs = new int[n], rs = new int[n];
		    for(int i = 0; i<n; i++)
		    {
			    gs[i] = Integer.parseInt(gFreq[i]);
			    rs[i] = Integer.parseInt(rFreq[i]);
		    }
		    int numSignals = 0;
		    boolean signal = false;
		    for(int i = 0; i<n; i++)
		    {
			    if(gs[i] == 0 && rs[i] >= t)
			    {
			        numSignals++;
			        signal = true;
			    }
		    }
		    double within = getWithinUniqueness(seq, seq.length() - n + 1);
		    if(freqMap == null)
		    {
		        freqMap = freqMap(allseqs, seq.length() - n + 1);
		    }
		    double external = getExternalUniqueness(freqMap, seq, seq.length() - n + 1);
		    if(call(signal, numSignals, n, within, external))
		    {
		        outs[p].println(seq);
		        outs[p].println(numSignals+" "+n);
		        outs[p].println(within);
		        outs[p].println(external);
		        count++;
		    }
	    }
	    offset += tot;
   	    System.out.println();
	    System.out.println(fns[p]);
	    System.out.println("Threshold for signal = " + t);
	    System.out.println("Number of SVs = " + tot);
	    System.out.println("Number with signal = " + count);

	    outs[p].close();
	}
}
static HashMap<String, Integer> freqMap(ArrayList<String> seqs, int k)
{
    HashMap<String, Integer> res = new HashMap<String, Integer>();
    for(String s : seqs)
    {
        for(int i = 0; i+k <= s.length(); i++)
        {
            String sub = s.substring(i, i+k);
            res.put(sub, res.containsKey(sub) ? (1 + res.get(sub)) : 1);
        }
    }
    return res;
}
static ArrayList<String> getAllSeqs(String fn) throws IOException
{
    ArrayList<String> res = new ArrayList<String>();
    Scanner input = new Scanner(new FileInputStream(new File(fn)));
    while(input.hasNext())
    {
        String line = input.nextLine();
   	    if(line.startsWith("Adding")) continue;
   	    if(line.charAt(0) < 'A' || line.charAt(0) > 'T') continue;
   	    res.add(line);
   	    for(int i = 0; i<3; i++) input.nextLine();
    }
    return res;
}
static boolean call(boolean signal, int numSignals, int n, double within, double external)
{
    return signal && within >= .5 && external >= .5;
}
static double getExternalUniqueness(HashMap<String, Integer> freqMap, String s, int k)
{
    int n = s.length();
    if(n < k || k <= 0)
    {
        return 1;
    }
    HashMap<String, Integer> map = new HashMap<String, Integer>();
    for(int i = 0; i+k <= s.length(); i++)
    {
        String sub = s.substring(i, i+k);
        map.put(sub, map.containsKey(sub) ? (1 + map.get(sub)) : 1);
    }
    int count = 0;
    for(String cur : map.keySet())
    {
        count += map.get(cur).equals(freqMap.get(cur)) ? map.get(cur) : 0;
    } 
    return 1.0 * count / (n - k + 1);
}
static double getWithinUniqueness(String s, int k)
{
    int n = s.length();
    if(n < k)
    {
        return 1;
    }
    HashSet<String> set = new HashSet<String>();
    for(int i = 0; i+k<=s.length(); i++)
    {
        set.add(s.substring(i, i+k));
    }
    return 1.0 * set.size() / (n - k + 1);
}
}
