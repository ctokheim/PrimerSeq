Java
====

These java files are the source code for the jar files located up one directory.
They are intended to add speed gains to processing the sequencing data compared
to python. They use the Sam-JDK from [Picard](http://picard.sourceforge.net/) and
the [BigWig API](http://code.google.com/p/bigwig/) from the [Broad Institute](http://www.broadinstitute.org/igv).

Files
=====

* Convert2SortedBam.java - Converts a SAM/BAM file to a sorted BAM file
* ExtractSamRegion.java - Extracts reads from a specified region of a BAM file
* ExtractBigRegion.java - Extract region from BigWig or BigBed file
