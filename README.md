PrimerSeq
=========

(c) 2012 Collin Tokheim GPLv3

http://primerseq.sf.net

Design RT-PCR primers using RNA-Seq data. 
Annotations often contain isoforms that only occur in 
certain tissues or cells. PrimerSeq tackles this issue by weighting 
isoforms supported in your RNA-Seq data when finding
flanking constitutive exons to place primers around a target.

You can find all of the licenses mentioned in this README in the `help/licenses` directory.

Installation
============

You will need [Java](http://www.oracle.com/technetwork/java/javase/downloads/java-se-jre-7-download-432155.html)
installed on your computer.

You will need [python 2.7.X](http://www.python.org/download/releases/2.7/) to run PrimerSeq.
If you are using Windows and are unfamiliar with python then I suggest you download the following
unofficial packages from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/):

* numpy
* matplotlib
* networkx
* wxPython

You will also need to install the latest version of [pygr](http://code.google.com/p/pygr/downloads/list).

Python Modules (Advanced)
=========================

If you are more familiar with python than you can install the python modules by using pip.  

```bash
$ pip install numpy
$ pip install matplotlib
$ pip install networkx
$ pip install pygr
```

Make sure networkx is version 1.7.0 otherwise PrimerSeq will not work.
If you wish to use the GUI you will also need wxPython.

```bash
$ pip install wxpython
```

wxPython may be a pain to install so you could additionally try using
`apt-get` if you use Ubuntu. 

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

Problems?
=========

If you find issues with PrimerSeq first consider looking at the [FAQ page](http://primerseq.sf.net/faq.html).
If that does not solve your problem, consider using GitHub's [issue submission](https://github.com/ctokheim/PrimerSeq/issues) since it will allow other users to see what has/is being fixed.
If needed you may send an email to me at primerseq@gmail.com.

