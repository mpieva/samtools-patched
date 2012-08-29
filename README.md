samtools-patched
================

A copy of samtools (from samtools.sf.net) with lots of patches:

 * Enhanced filter syntax for 'samtools view'.  The '-f' and 'F'
   switches not take the same mnemonic codes as those used in the
   output.  In addition, the XF tagged field contains more flags, we use
   1(m) for merged and 2(t) for adapter-trimmed reads, Filtering on the
   sequence length is also possible ('-m' or '-M').

 * More statistics in 'samtools flagstat'.  Also, the same statistics are
   available in table form per read-group through 'samtools flagstatx'.

 * The new command 'samtools covstat' tallies coverage of the target by
   read-group and by target sequence.
   
 * 'samtools view' can produce FastQ as output format.  Alignments are
   not included, but most common conventions are followed.  (This code
   avoids producing special characters at the beginning of quality line
   to make life easier for weak FastQ parsers.  Please don't use FastQ
   if you have the choice.)

 * Cleaner handling of read groups in commands that respect them.
   Single characters are generally allowed as valid read groups.

 * More error checking and reporting, especially on out-of-disk-space
   conditions.

