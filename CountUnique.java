/*
 * Given a set of insertions, a reference, and short-read data, counts how many times 
 * each kmer in the insertions occurs in both the reference and set of reads
 */
import java.util.*;
import java.io.*;
public class CountUnique {
	static int k = 30;
public static void main(String[] args) throws IOException
{
	String fastaFn = args[0];
	String vcfFn = args[1];
	String read1Fn = args[2];
	String read2Fn = args[3];
	Scanner vcfInput = new Scanner(new FileInputStream(new File(vcfFn)));
	Scanner fastaInput = new Scanner(new FileInputStream(new File(fastaFn)));
	Scanner read1Input = new Scanner(new FileInputStream(new File(read1Fn)));
	Scanner read2Input = new Scanner(new FileInputStream(new File(read2Fn)));
	ArrayList<String> seqs = new ArrayList<String>();
	while(vcfInput.hasNext())
	{
		String line = vcfInput.nextLine();
		if(line.startsWith("#")) continue;
		if(!line.contains("SEQ=")) continue;
		String sub = line.substring(line.indexOf("SEQ=")+4);
		sub = sub.substring(0,sub.indexOf(';'));
		seqs.add(sub);
	}
	long[][] g = new long[seqs.size()][];
	for(int i = 0; i<seqs.size(); i++)
	{
		if(seqs.get(i).length() < k)
		{
			g[i] = new long[0];
			continue;
		}
		g[i] = new long[seqs.get(i).length() - k + 1];
		ArrayList<Long> curKmers = kmerize(seqs.get(i));
		for(int j = 0; j<g[i].length; j++) g[i][j] = curKmers.get(j);
	}
	HashMap<Long, ArrayList<Integer>> map = new HashMap<Long, ArrayList<Integer>>();
	for(int i = 0; i<g.length; i++)
		for(int j = 0; j<g[i].length; j++)
		{
			long cur = g[i][j];
			if(!map.containsKey(cur)) map.put(cur, new ArrayList<Integer>());
			map.get(cur).add(j * g.length + i);
		}
	int[][] counts = new int[g.length][];
	for(int i = 0; i<g.length; i++) counts[i] = new int[g[i].length];
	String name = fastaInput.nextLine();
	StringBuilder sb = new StringBuilder("");
	while(fastaInput.hasNext())
	{
		String cur = fastaInput.nextLine();
		if(cur.charAt(0) == '>')
		{
			System.out.println("Adding " + name);
			process(counts, map, sb.toString());
			sb = new StringBuilder("");
			name = cur;
		}
		else sb.append(cur);
	}
	process(counts, map, sb.toString());
	
	int[][] counts2 = new int[g.length][];
	for(int i = 0; i<g.length; i++) counts2[i] = new int[g[i].length];
	while(read1Input.hasNext())
	{
		read1Input.nextLine();
		String a = read1Input.nextLine();
		read2Input.nextLine();
		String b = read2Input.nextLine();
		for(int i = 0; i<2; i++)
		{
			read1Input.nextLine();
			read2Input.nextLine();
		}
		process(counts2, map, a);
		process(counts2, map, b);
	}
	int count = 0;
	for(int i = 0; i<seqs.size(); i++)
	{
		String s = seqs.get(i);
		System.out.println(s);
		int min = 987654321;
		for(int j = 0; j<g[i].length; j++)
		{
			min = Math.min(min, counts[i][j]);
			System.out.print(counts[i][j]+((j == g[i].length - 1) ? "\n" : " "));
		}
		for(int j = 0; j<g[i].length; j++)
		{
			System.out.print(counts2[i][j]+((j == g[i].length - 1) ? "\n" : " "));
		}
		System.out.println(min);
	}
	System.out.println(count);
}
static void process(int[][] counts, HashMap<Long, ArrayList<Integer>> map, String s)
{
	long hash = 0;
	int[] vals = new int[256];
	vals['A'] = 0;
	vals['C'] = 1;
	vals['G'] = 2;
	vals['T'] = 3;
	for(int i = 0; i<k; i++) hash = (hash << 2) + vals[s.charAt(i)];
	if(map.containsKey(hash))
	{
		for(int x : map.get(hash))
		{
			counts[x%counts.length][x/counts.length]++;
		}
	}
	for(int i = k; i<s.length(); i++)
	{
		hash = hash & ((1L << (2 * (k - 1))) - 1);
		hash = (hash << 2) + vals[s.charAt(i)];
		if(map.containsKey(hash))
		{
			for(int x : map.get(hash))
			{
				counts[x%counts.length][x/counts.length]++;
			}
		}
	}
}
static ArrayList<Long> kmerize(String s)
{
	ArrayList<Long> res = new ArrayList<Long>();
	long hash = 0;
	int[] vals = new int[256];
	vals['A'] = 0;
	vals['C'] = 1;
	vals['G'] = 2;
	vals['T'] = 3;
	for(int i = 0; i<k; i++) hash = (hash << 2) + vals[s.charAt(i)];
	res.add(hash);
	for(int i = k; i<s.length(); i++)
	{
		hash = hash & ((1L << (2 * (k - 1))) - 1);
		hash = (hash << 2) + vals[s.charAt(i)];
		res.add(hash);
	}
	return res;
}
static void add(HashMap<Long, Integer> kmers, String s)
{
	if(s.length() < k) return;
	long hash = 0;
	int[] vals = new int[256];
	vals['A'] = 0;
	vals['C'] = 1;
	vals['G'] = 2;
	vals['T'] = 3;
	for(int i = 0; i<k; i++) hash = (hash << 2) + vals[s.charAt(i)];
	kmers.put(hash, kmers.containsKey(hash) ? (1 + kmers.get(hash)) : 1);
	for(int i = k; i<s.length(); i++)
	{
		hash = hash & ((1L << (2 * (k - 1))) - 1);
		hash = (hash << 2) + vals[s.charAt(i)];
		kmers.put(hash, kmers.containsKey(hash) ? (1 + kmers.get(hash)) : 1);
	}
}
}
