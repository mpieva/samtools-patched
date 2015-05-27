CC=			gcc
# CFLAGS=		-ggdb -Wall -O0 #-m64 #-arch ppc
CFLAGS=		-Wall -O2 -ggdb #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=1 \
			-DGIT_VERSION="$(shell git describe --always)"
KNETFILE_O=	knetfile.o
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o bam_md.o razf.o faidx.o bedidx.o \
			$(KNETFILE_O) bam_sort.o sam_header.o bam_reheader.o kprobaln.o bam_cat.o
AOBJS=		bam_tview.o bam_plcmd.o sam_view.o	\
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o	\
			bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o \
			cut_target.o phase.o bam2depth.o
PROG=		samtools razip faidx bgzip
INCLUDES=	-I.
SUBDIRS=	. bcftools misc
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

install: all-recur
	install -d -m 755 $(prefix)/bin
	install -d -m 755 $(prefix)/include/bam
	install -d -m 755 $(prefix)/lib
	install -d -m 755 $(prefix)/share/man/man1
	-install -s -m 755 -b samtools $(prefix)/bin
	-install -s -m 755 -b samtools $(prefix)/bin/sam
	install -s -m 755 -b razip $(prefix)/bin
	install -s -m 755 -b bgzip $(prefix)/bin
	install -s -m 755 -b faidx $(prefix)/bin
	install    -m 644 -b libbam.a $(prefix)/lib
	install    -m 644    *.h $(prefix)/include/bam/
	-install -s -m 755 -b bcftools/bcftools $(prefix)/bin
	-install    -m 755    bcftools/vcfutils.pl $(prefix)/bin
	install    -m 644    samtools.1 $(prefix)/share/man/man1
	-install -s -m 755 -b misc/seqtk $(prefix)/bin
	-install -s -m 755 -b misc/wgsim $(prefix)/bin



.PHONY:all lib clean cleanlocal install
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

samtools:lib-recur $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) -Lbcftools $(LIBPATH) libbam.a -lbcf $(LIBCURSES) -lm -lz
		cp samtools sam

razip:razip.o razf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ razf.o razip.o $(KNETFILE_O) -lz

bgzip:bgzip.o bgzf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ bgzf.o bgzip.o $(KNETFILE_O) -lz

faidx: faidx.c razf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ -DFAIDX_MAIN faidx.c razf.o $(KNETFILE_O) -lz

libbam.1.dylib-local:$(LOBJS)
		libtool -dynamic $(LOBJS) -o libbam.1.dylib -lc -lz

libbam.so.1-local:$(LOBJS)
		$(CC) -shared -Wl,-soname,libbam.so -o libbam.so.1 $(LOBJS) -lc -lz

dylib:
		@$(MAKE) cleanlocal; \
		case `uname` in \
			Linux) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.so.1-local;; \
			Darwin) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.1.dylib-local;; \
			*) echo 'Unknown OS';; \
		esac


cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib

clean:cleanlocal-recur

bam.o: bam.c bam.h bgzf.h bam_endian.h kstring.h sam_header.h
bam2bcf.o: bam2bcf.c bam.h bgzf.h kstring.h bam2bcf.h errmod.h \
 bcftools/bcf.h
bam2bcf_indel.o: bam2bcf_indel.c bam.h bgzf.h bam2bcf.h errmod.h \
 bcftools/bcf.h kaln.h kprobaln.h khash.h ksort.h
bam2depth.o: bam2depth.c bam.h bgzf.h
bam_aux.o: bam_aux.c bam.h bgzf.h khash.h
bam_cat.o: bam_cat.c bgzf.h bam.h
bam_color.o: bam_color.c bam.h bgzf.h
bam_import.o: bam_import.c kstring.h bam.h bgzf.h sam_header.h kseq.h \
 khash.h
bam_index.o: bam_index.c khash.h ksort.h bam_endian.h bam_index.h bam.h \
 bgzf.h
bam_lpileup.o: bam_lpileup.c bam.h bgzf.h ksort.h
bam_mate.o: bam_mate.c bam.h bgzf.h
bam_md.o: bam_md.c faidx.h sam.h bam.h bgzf.h kstring.h kaln.h kprobaln.h
bam_pileup.o: bam_pileup.c sam.h bam.h bgzf.h
bam_plcmd.o: bam_plcmd.c sam.h bam.h bgzf.h faidx.h kstring.h bam2bcf.h \
 errmod.h bcftools/bcf.h sample.h
bam_reheader.o: bam_reheader.c bgzf.h bam.h
bam_rmdup.o: bam_rmdup.c sam.h bam.h bgzf.h khash.h
bam_rmdupse.o: bam_rmdupse.c sam.h bam.h bgzf.h khash.h klist.h
bam_sort.o: bam_sort.c ksort.h bam_stat.h bam.h bgzf.h khash.h \
 bam_index.h
bam_stat.o: bam_stat.c bam_stat.h bam.h bgzf.h khash.h
bam_tview.o: bam_tview.c
bamtk.o: bamtk.c bam.h bgzf.h
bedidx.o: bedidx.c ksort.h kseq.h khash.h
bgzf.o: bgzf.c bgzf.h khash.h
bgzip.o: bgzip.c bgzf.h
cut_target.o: cut_target.c bam.h bgzf.h errmod.h faidx.h
errmod.o: errmod.c errmod.h ksort.h
faidx.o: faidx.c faidx.h khash.h razf.h
kaln.o: kaln.c kaln.h
knetfile.o: knetfile.c knetfile.h
kprobaln.o: kprobaln.c kprobaln.h
kstring.o: kstring.c kstring.h
phase.o: phase.c bam.h bgzf.h errmod.h kseq.h khash.h ksort.h
razf.o: razf.c razf.h
razip.o: razip.c razf.h
sam.o: sam.c faidx.h sam.h bam.h bgzf.h
sam_header.o: sam_header.c sam_header.h khash.h
sam_view.o: sam_view.c sam_header.h sam.h bam.h bgzf.h faidx.h kstring.h \
 khash.h
sample.o: sample.c sample.h kstring.h khash.h
