import java.util.*;
import java.io.*;
public class VCFDiff {
	static int INTERSECT = 1;
	static int DIFF = 0;
	static int setting = DIFF;
	static boolean filter = true; // Whether or not to take only lines with SVTYPE=INS
	static int maxDist = 10;
	static double similarity = 0;
	static int INSERT = 0;
	static int DELETE = 1;
	static int svType = INSERT;
	static double lengthThreshold = 0.5;
public static void main(String[] args) throws IOException
{
	String fn1 =args[0];
	String fn2 = args[1];
	String outFn = args[2];
	setting = args[3].equals("intersect") ? INTERSECT : DIFF;
	svType = args[4].equals("insert") ? INSERT : DELETE;
	maxDist = Integer.parseInt(args[5]);
	similarity = Double.parseDouble(args[6]);
	diff(fn1, fn2, outFn);
}
static void diff(String fn1, String fn2, String outFn) throws IOException
{
	TreeSet<Insertion> sv1 = parseFile(fn1);
	TreeSet<Insertion> sv2 = parseFile(fn2);
	TreeSet<Insertion> diff = diffInsertions(sv1, sv2);
	PrintWriter out = new PrintWriter(new File(outFn));
	printHeader(fn1, out);
	filterVCF(fn1, diff, out);
	out.close();
}
/*
 * Prints all SVs in a VCF file which are also present in a given set
 */
static void filterVCF(String fn, TreeSet<Insertion> svs, PrintWriter out) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.charAt(0) == '#') continue;
		if(filter && !line.contains((svType == DELETE) ? "SVTYPE=DEL" : "SVTYPE=INS")) continue;
		Insertion cur = new Insertion(line);
		if(svs.contains(cur)) out.println(line);
	}
}
/*
 * Prints the header of a VCF file
 */
static void printHeader(String fn, PrintWriter out) throws IOException
{
	@SuppressWarnings("resource")
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0) continue;
		if(line.charAt(0) == '#') out.println(line);
		else break;
	}
}
/*
 * Returns SV's in sv1 but not in sv2
 */
static TreeSet<Insertion> diffInsertions(TreeSet<Insertion> sv1, TreeSet<Insertion> sv2)
{
	TreeSet<Insertion> res = new TreeSet<Insertion>();
	for(Insertion cur : sv1)
	{
		// Check nearby insertions in sv2 for sufficient similarity
		boolean found = false;
		Insertion ceil = sv2.ceiling(cur);
		while(!found && ceil != null && cur.chr.equals(ceil.chr) && Math.abs(cur.pos - ceil.pos) <= maxDist)
		{
			if(cur.similar(ceil))
			{
				found = true;
				break;
			}
			ceil = sv2.higher(ceil);
		}
		Insertion floor = sv2.floor(cur);
		while(!found && floor != null && cur.chr.equals(floor.chr) && Math.abs(cur.pos - floor.pos) <= maxDist)
		{
			if(cur.similar(floor))
			{
				found = true;
				break;
			}
			floor = sv2.lower(floor);
		}
		if((!found && setting == DIFF) || (found && setting == INTERSECT))
	    {   
	        res.add(cur);
	    }
	}
	return res;
}
/*
 * Reads and indexes a list of insertions from a VCF file
 */
static TreeSet<Insertion> parseFile(String fn) throws IOException
{
	TreeSet<Insertion> res = new TreeSet<>();
	@SuppressWarnings("resource")
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.charAt(0) == '#') continue;
		if(filter && svType == DELETE && !line.contains("SVTYPE=DEL")) continue;
		if(filter && svType == INSERT && !line.contains("SVTYPE=INS")) continue;
		Insertion cur = new Insertion(line);
		if(cur.length == 0) continue;
		res.add(new Insertion(line));
	}
	return res;
}
/*
 * Extracts the sequence of an insertion from its VCF entry
 */
static String getSeq(String line)
{
	String pattern = "SEQ=";
	int idx = line.indexOf(pattern);
	if(idx == -1) return "";
	line = line.toUpperCase();
	idx += pattern.length();
	int end = idx;
	while(end < line.length() && isBasePair(line.charAt(end))) end++;
	return line.substring(idx, end);
}
/*
 * Extracts the sv length of an insertion from its VCF entry
 */
static int getLength(String line)
{
	String pattern = "SVLEN=";
	int idx = line.indexOf(pattern);
	if(idx == -1) return 0;
	line = line.toUpperCase();
	idx += pattern.length();
	int end = idx;
	int length = 0;
	while(end < line.length() && line.charAt(end) >= '0' && line.charAt(end) <= '9')
	{
	    length = length * 10 + (line.charAt(end) - '0');
	    end++;
	}
	return length;
}
/*
 * Returns whether or not a character is a base pair
 */
static boolean isBasePair(char c)
{
	return c >= 'A' && c <= 'Z';
}
/*
 * Returns longest common subsequence of two strings s and t
 */
static int lcs(String s, String t)
{
	int n = s.length(), m = t.length();
	int[][] dp = new int[n+1][m+1];
	for(int i = 1; i<=n; i++)
		for(int j = 1; j<=m; j++)
		{
			dp[i][j] = Math.max(dp[i-1][j], dp[i][j-1]);
			int cur = (s.charAt(i-1) == t.charAt(j-1)) ? 1 : 0;
			dp[i][j] = Math.max(dp[i][j], cur + dp[i-1][j-1]);
		}
	return dp[n][m];
}
/*
 * Represents an insertion by its position and sequence
 */
static class Insertion implements Comparable<Insertion>
{
	String chr;
	int pos;
	String seq;
	int length;
	Insertion(String line)
	{
		String[] tokens = line.split("\t");
		chr = tokens[0];
		if(chr.startsWith("chr")) chr = chr.substring(3);
		if(chr.equals("MT")) chr = "M";
		pos = Integer.parseInt(tokens[1]);
		seq = getSeq(line);
		length = getLength(line);
	}
	public int compareTo(Insertion o)
	{
		if(!chr.equals(o.chr)) return chr.compareTo(o.chr);
		if(pos != o.pos) return pos - o.pos;
		return seq.compareTo(o.seq);
	}
	boolean similar(Insertion o)
	{
	    if(svType == INSERT)
	    {
		    String s = seq, t = o.seq;
		    int minLength = Math.min(s.length(), t.length());
		    int maxLength = Math.max(s.length(), t.length());
		    if(minLength < similarity * maxLength) return false;
		    return lcs(s, t) >= similarity * minLength;
		}
		else
		{
		    if(length < lengthThreshold * o.length) return false;
		    if(o.length < length * lengthThreshold) return false;
		    return true;
		}
	}
}
}
