import java.util.*;
import java.io.*;

// This class grabs the gene id and tx id from the attribute field in a GTF
class AttributeField {
	String gene_id;
	String tx_id;
	String text;
	
	public AttributeField(String attr){
		this.text = attr.trim();  // store input for toString method
		String[] split_string = this.text.split("[ ]*;*[ ]+");
		int i; // counter
		String tmp, next_str;
		for(i=0;i<split_string.length-1;i+=2){
			tmp = split_string[i];  // shorter variable name
			next_str = split_string[i+1];  // the value of the associated key
			if (tmp.toLowerCase().equals("gene_id")){
				this.gene_id = next_str.replaceAll("\"", "");
			}else if (tmp.toLowerCase().equals("transcript_id")){
				this.tx_id = next_str.replaceAll("\"", "");
			}
		}
	}
	
	public String getGeneID(){ return this.gene_id; }
	
	public String getTxID(){ return this.tx_id; }
	
	// output the unmodified text
	public String toString(){
		return this.text;
	}
}

// Puts a String array and a AttributeField object into one single object for comparator
// Perhaps not the best you of a class but it avoids type casting issues when inheriting from comparator
class LineData {
	String[] line;
	AttributeField attr;
	
	public LineData(String[] splitted_line) {
		this.line = splitted_line;
		this.attr = new AttributeField(this.line[8]);
		this.line[8] = this.attr.toString();
	}
	
	// getter for string array
	public String[] getLine(){
		return this.line;
	}
	
	// getter for attr
	public AttributeField getAttr(){
		return this.attr;
	}
}

class SortGtfComparator implements Comparator<Object> {
	
	public int compare(Object o1, Object o2){
		LineData l1 = (LineData) o1;
		LineData l2 = (LineData) o2;
		String[] line1 = l1.getLine();
		String[] line2 = l2.getLine();
		
		// compare chromosome names
		String seqname1 = (String) line1[0];
		int seqCompare = seqname1.compareTo((String) line2[0]);
		
		// compare gene/tx ids
		AttributeField attr1 = l1.getAttr();
		AttributeField attr2 = l2.getAttr();
		try{
			int geneCompare = attr1.getGeneID().compareTo(attr2.getGeneID());
			int txCompare = attr1.getTxID().compareTo(attr2.getTxID());
			
			// comparison logic
			if(seqCompare < 0){
				return -1;
			} else if (seqCompare > 0){
				return 1;
			} else if (geneCompare < 0) {
				return -1;
			} else if (geneCompare > 0) {
				return 1;
			} else if (txCompare < 0){
				return -1;
			} else if (txCompare > 0){
				return 1;
			} else if (Integer.parseInt(line1[3]) < Integer.parseInt(line2[3])){
				return -1;
			} else if (Integer.parseInt(line1[3]) > Integer.parseInt(line2[3])){
				return 1;
			} else if (Integer.parseInt(line1[4]) < Integer.parseInt(line2[4])){
				return -1;
			} else if (Integer.parseInt(line1[4]) > Integer.parseInt(line2[4])){
				return 1;
			} else {
				return 0;
			}
		} catch(Exception e){
			System.out.println(Arrays.toString(line1));
			System.out.println(Arrays.toString(line2));
			e.printStackTrace();
			System.exit(1);
			return 0;
		}
	}
}




// Read in tab delimited gtf file
class TabDelimReader {
	// ArrayList<String[]> lines = new ArrayList<String[]>();
	ArrayList<LineData> lines = new ArrayList<LineData>();
	public TabDelimReader(String fname){
		try {
			String[] tmp_line;
			Scanner scan = new Scanner(new File(fname));
			while(scan.hasNextLine()){
				tmp_line = scan.nextLine().split("\t");  // hold line data temporarily
				
				// Ignore lines that are not "exon" features for a speed up in performance
				if(tmp_line[2].equals("exon")) {
					this.lines.add(new LineData(tmp_line));
				}
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public ArrayList<LineData> getLines(){
		return this.lines;
	}
}


// SortGtf is the main class that drives sorting of a GTF by [seqname, gene id, tx id, start, end]
public class SortGtf {

	// The join method acts just like python's string.join()
	public static String join(String[] fields, String separator){
	    StringBuilder stringBuilder = new StringBuilder();
	    for(int i = 0; i < fields.length; i++){
            stringBuilder.append(fields[i]);
            if(i < fields.length - 1){
            	stringBuilder.append(separator);
            }
	    }
	    return stringBuilder.toString();                           
	}
	
	
	public static void main(String[] args){
		// complain about args
		if (args.length != 2){
			System.out.println("SortGtf.jar requires exactly two arguments.");
			System.out.println("Usage: java -jar SortGtf.jar input.gtf output.gtf");
			System.exit(1);
		}
		
		// Read in and then sort the GTF
		TabDelimReader inputReader = new TabDelimReader(args[0]);
		ArrayList<LineData> lines = inputReader.getLines();  // read in GTF (tab-separated)
		Object[] data = lines.toArray();
		Arrays.sort(data, new SortGtfComparator()); // sort the data
		
		// output the sorted gtf
		try{
			Formatter output=new Formatter(args[1]);
			
			String[] dataLine;
			int i;
			LineData tmp;
			for(i=0;i<data.length;i++){
				// dataLine = (String[]) data[i];
				tmp = (LineData) data[i];
				dataLine = (String[]) tmp.getLine();
				output.format("%s\n", SortGtf.join(dataLine, "\t"));
			}
			output.close();
		} catch(FileNotFoundException e){
			e.printStackTrace();
		}
		
	}
}
