hapflk_pool: FLK test for Pool-Seq data
=======================================

## Motivation

__FLK__ ([Bonhomme _et al._ 2010](http://www.genetics.org/content/186/1/241.abstract)) is a test for selection scans based on population allele frequencies.
Nevertheless, the current implementation in python ([hapFLK](https://forge-dga.jouy.inra.fr/projects/hapflk)) was develop for using haplotypes (due to the
extension of the method by [Fariello _et al._ 2013](http://www.genetics.org/content/193/3/929.abstract)) in the PLINK format.

A new feature implemented in hapFLK release 1.3.0 is the tree rooting using a MLE instead where an outgroup is not available (more details 
[here](https://forge-dga.jouy.inra.fr/projects/hapflk/news)). This is extremely important if the test is performed in two populations, and also to avoid the 
bias in the kinship matrix as shown in [Gautier 2015](http://www.genetics.org/content/201/4/1555.long).

Although there is an [R package](https://qgsp.jouy.inra.fr/index.php?option=com_content&view=article&id=50&Itemid=55) for compute FLK directly from allele frequencies,
it lacks the option to use the MLE for rooting. 

Thus, I changed the source code from hapFLK realease 1.3 to include the option to use the software with a widely used format in Pool-Seq data: synchronized population
files (sync files) as described in [PoPoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/), including gzipped versions. Like that, the software here allows 
to compute FLK for Pool-Seq data with the hapFLK software including the latest tree rooting algorithm, without removing options from use the haplotype model for other 
datasets.

## hapFLK information

For more details about hapFLK check the [README.txt](https://github.com/magicDGS/hapFLK-extension/blob/master/README.txt) from the original code and previous references.