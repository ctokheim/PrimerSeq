PrimerSeq
=========

(c) 2012-2013 Collin Tokheim GPLv3

http://primerseq.sf.net

Design RT-PCR primers that evaluate the **P**ercent **S**pliced **I**n (PSI) metric using RNA-Seq data.
Thus differential Alternative Splicing (AS) events can quickly be validated from programs like [MATS](http://rnaseq-mats.sourceforge.net/).
Annotations often contain isoforms that only occur in
certain tissues or cells. PrimerSeq tackles this issue by weighting
isoforms supported in your RNA-Seq data when finding
flanking constitutive exons to place primers around a target.

You can find all of the licenses mentioned in this README in the `help/licenses` directory.

Installation
============

You will need [Java](http://www.oracle.com/technetwork/java/javase/downloads/java-se-jre-7-download-432155.html)
installed on your computer regardless if you are on Windows or Linux.

You can download an installer for Windows [here](http://sourceforge.net/projects/primerseq/files/PrimerSeq/). There is a step-by-step [installation guide](http://primerseq.sourceforge.net/windows.html) for Windows. Installation for linux will need to be from source. To install from source on Ubuntu (>=12.04) or Linux Mint 13 follow the guide [here](http://primerseq.sourceforge.net/linux.html).

Installation from Source
========================

If you are using Ubuntu (>=12.04) or Linux Mint 13, please follow the [linux installation guide](http://primerseq.sourceforge.net/linux.html) otherwise you will need to perform installation your self.

You will need [python 2.7.X](http://www.python.org/download/releases/2.7/) to run PrimerSeq. For linux
users, you can install python 2.7.3 (also with numpy/scipy) in your ```~/local/bin``` folder by using this [script](https://gist.github.com/4507404).
If you are using Windows and are unfamiliar with python then I suggest you download the following
unofficial packages from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/):

* numpy
* matplotlib
* networkx
* wxPython

You will also need to install the latest version of [pygr](http://code.google.com/p/pygr/downloads/list).

Python Modules (Advanced)
=========================

View the [install.ubuntu.sh](https://github.com/ctokheim/PrimerSeq/blob/master/install.ubuntu.sh) script for actual commands necessary to install PrimerSeq on ubuntu. The *install.ubuntu.sh* script should install everything required except Java.

Make sure networkx is version 1.7.0 otherwise PrimerSeq will not work.

Primer3
=======

Primer3 does the actual primer design after appropriate sequences are extracted.
Primer3 is assumed to be in the `PrimerSeq/primer3`
directory. That means primer3_core should be found `PrimerSeq/primer3/src/primer3_core`.
A custom path for Primer3 can be specified by modifying PrimerSeq.cfg. Primer3 is licensed
under GPLv2.0.

Java
====

The Java JRE should also be installed since jar files in the bin directory
are used to handle the sequencing data.

The Jar files include external libraries which are included with out any modification.
The external libraries are listed below with their license version:

* [Sam-JDK](http://picard.sourceforge.net/) handles BAM/SAM files ([Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0.html), [MIT](http://opensource.org/licenses/MIT))
* [BigWig API](http://code.google.com/p/bigwig/) extracts regions from a BigWig file ([LGPL v2.1](http://www.gnu.org/licenses/lgpl-2.1.html))

If you only intend to run PrimerSeq then you do **NOT** need to obtain the external Java APIs since they are bundled with the jar files found in the bin directory.

Problems?
=========

If you find issues with PrimerSeq first consider looking at the [FAQ page](http://primerseq.sf.net/faq.html).
If that does not solve your problem, consider using GitHub's [issue submission](https://github.com/ctokheim/PrimerSeq/issues) since it will allow other users to see what has/is being fixed.
If needed you may send an email to me at primerseq@gmail.com.

