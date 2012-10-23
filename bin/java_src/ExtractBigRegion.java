import java.io.*;
import java.util.*;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BBFileHeader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;


public class ExtractBigRegion {
	
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


	public static void main(String[] args) throws IOException{
		if(args.length != 6){
			System.out.println("This utility extracts a region from a BigWig or BigBed file into their corresponding wig or bed file.\n");
			System.out.println("Usage:\njava -jar ExtractBigRegion.jar <input.(bigWig|bigBed)> <output.(wig|bed)> <chr> <start> <end> <contained>\n");
			System.out.println("Example:\njava -jar ExtractBigRegion.jar example.bigWig out.wig chr1 2000 30000 true\n");
			System.out.println("NOTE: The contained parameter is a boolean (true or false). \"true\" means features with only partial overlap are not extracted.\n");
			System.exit(0);
		}
		
		// parse input
		int start = Integer.parseInt(args[3]);
		int end = Integer.parseInt(args[4]);
		if(end < start){
			// flip start and stop if necessary
			int tmp = start;
			start = end;
		    end = tmp;
		}
		boolean contained = false;
		if(args[5] == "true"){
			contained = true; //set to true
		}
		
		BBFileReader reader=new BBFileReader(args[0]); // get big file reader
		BBFileHeader header = reader.getBBFileHeader();  // get big file header
		
		// check header
		if(!header.isHeaderOK()){
			System.out.println("The header is not OK!");
			System.exit(1);
		}
		
		// make sure header and file agree on file type
		if(!(reader.isBigBedFile() == header.isBigBed() && reader.isBigWigFile() == header.isBigWig())){
			System.out.println("The file and headers do not agree on file type");
			System.exit(1);
		}
		
		// Complain if neither BigBed nor BigWig
		//if(!(header.isBigBed() || header.isBigWig())){
		//	System.out.println("Yikes! This only works for BigBed and BigWig");
		//	System.exit(1);
		//}
		
		Formatter output=new Formatter(args[1]);  // file to write output
		
		// output information according to file type
		if(header.isBigWig()){
			BigWigIterator iter = reader.getBigWigIterator(args[2], start, args[2], end, contained);
			while(iter.hasNext()){
				WigItem wigItem = iter.next();
				output.format("%s\t%d\t%d\t%d\n", wigItem.getChromosome(), wigItem.getStartBase(), wigItem.getEndBase(), (int) wigItem.getWigValue());
			}
		}else if(header.isBigBed()){
			BigBedIterator iter = reader.getBigBedIterator(args[2], start, args[2], end, contained);
			while(iter.hasNext()){
				BedFeature bedFeature = iter.next();
				output.format("%s\t%d\t%d\t%s\n", bedFeature.getChromosome(), bedFeature.getStartBase(), bedFeature.getEndBase(),
							  ExtractBigRegion.join(bedFeature.getRestOfFields(), "\t"));  // the join method is similar to pythons '<sep>'.join
			}
		}

		output.close();  // make sure to close file
	}
}
