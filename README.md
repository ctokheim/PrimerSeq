PrimerSeq
=========

Design RT-PCR primers using RNA-Seq data. 
Annotations often contain isoforms that only occur in 
certain tissues or cells. PrimerSeq tackles this issue by weighting 
isoforms supported in your RNA-Seq data when finding
flanking constitutive exons to place primers around a target.

Note
====

Primer3 is required and is assumed to be in the `PrimerSeq/primer3`
directory. That means primer3_core should be found `PrimerSeq/primer3/src/primer3_core`.
The Java JRE should also be installed since jar files in the bin directory
are used to handle the sequencing data.
