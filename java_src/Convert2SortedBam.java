import net.sf.samtools.*; 
import java.io.File;

public class Convert2SortedBam {
	public static void main(String[] args){
		// explain command line args if they don't enter the correct number
		if(args.length != 2){
			System.out.println("This utility converts SAM/BAM file to sorted BAM file.\n");
			System.out.println("Usage:\njava -jar Convert2SortedBam.jar <input.(sam|bam)> <output.sorted.bam>\n");
			System.out.println("Example:\njava -jar Convert2SortedBam.jar example.sam example.sorted.bam");
			System.exit(0);
		}
		
		// check files
		if(!args[0].endsWith(".sam") && !args[0].endsWith(".bam")){
			System.out.println("The input file must be a .sam or .bam file");
			System.exit(0);
		}
		if(!args[1].endsWith(".bam")){
			System.out.println("The output must be a .bam file");
			System.exit(0);
		}
		
		// handle files
		File input = new File(args[0]);  // input sam file object
        File outputBamFile = new File(args[1]);
        SAMFileReader sam = new SAMFileReader(input);
        sam.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);  // specify that it should sort file
        final SAMFileWriter outputBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(sam.getFileHeader(), false, outputBamFile); // output
        
        // iterate over all reads
        SAMRecordIterator recIterator = sam.iterator(); // iterator over all SAMRecord's
        //recIterator.assertSorted(SAMFileHeader.SortOrder.coordinate);  // complain about not being sorted
        while(recIterator.hasNext()){
        	outputBam.addAlignment(recIterator.next());
        }
        
        // close files
        outputBam.close();
        sam.close();

	}
}
