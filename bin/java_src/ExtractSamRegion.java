import net.sf.samtools.*; 
import java.io.File;

public class ExtractSamRegion {
	public static void main(String[] args){
		// explain command line args if they don't enter the correct number
		if(args.length != 5){
			System.out.println("This utility gets the reads in a bam file that overlap with the 1-based coordinates specified by the user.\n");
			System.out.println("Usage:\njava -jar ExtractSamRegion.jar <input.bam> <output.sam> <chr> <start> <end>\n");
			System.out.println("Example:\njava -jar ExtractSamRegion.jar example.bam out.sam chr1 2000 30000");
			System.exit(0);
		}
		 
		// handle files
		if(!args[0].endsWith(".bam")){
			System.out.println("The input file must be a .bam file");
			System.exit(0);
		}
		File inputBam = new File(args[0]);
		if(!args[1].endsWith(".sam")){
			System.out.println("ExtractSamRegion outputs the region as a .sam file");
			System.exit(0);
		}
		File outputSamFile = new File(args[1]);
        SAMFileReader checkBam = new SAMFileReader(inputBam);
		
        // Make sure bam is actually binary and thus not sam
        if (!checkBam.isBinary()) {
        	System.out.println("The input file must be bam, not sam. The input file was not binary as expected.");
        	System.exit(1);
        }
        
        // index bam file if necessary
        if(!checkBam.hasIndex()){
            SAMFileReader tmpBam = new SAMFileReader(inputBam);
        	File indexFile = new File(args[0] + ".bai");
        	BAMIndexer index = new BAMIndexer(indexFile, tmpBam.getFileHeader());
        	tmpBam.enableFileSource(true); 
        	SAMRecordIterator myIterator = tmpBam.iterator();
        	SAMRecord t;
        	while(myIterator.hasNext()){
        		t = myIterator.next();
        		index.processAlignment(t);
        	}
        	index.finish();
        	tmpBam.close();
        }
        
        checkBam.close();  // close the bam file that was used to check if it was valid
        SAMFileReader bam = new SAMFileReader(inputBam);  // Read in bam after ensuring index
        
        // parse the start, end argument inputs
        int start = Integer.parseInt(args[3]);
        int end = Integer.parseInt(args[4]);
        if(start > end){
        	System.out.println("Yikes your start is greater than your end!");
        	System.exit(1);
        }
        
        
        final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(bam.getFileHeader(),
                true, outputSamFile); // output
        
        // iterate over reads in the specified region
        int totalReads = 0;
        SAMRecordIterator regionIterator = bam.queryOverlapping(args[2], start, end);
        SAMRecord tmp;
        while(regionIterator.hasNext()){
        	tmp = regionIterator.next();
        	outputSam.addAlignment(tmp);
        	totalReads++;
        }

        System.out.println("A total of " + Integer.toString(totalReads) + " were found overlapping " + args[2] + ":" + args[3] + "-" + args[4]);
        
        // close files
        outputSam.close();
        bam.close();

	}
}
