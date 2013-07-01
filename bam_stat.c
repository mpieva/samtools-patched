#include <unistd.h>
#include <assert.h>
#include "bam_stat.h"

#define flagstat_loop(s, c) do {										\
		int w = ((c)->flag & BAM_FQCFAIL)? 1 : 0;	 					\
		++(s)->n_reads[w];                                              \
		if ((c)->flag & BAM_FPAIRED) {                                  \
			++(s)->n_pair_all[w];                                       \
			if ((c)->flag & BAM_FPROPER_PAIR) ++(s)->n_pair_good[w];    \
			if ((c)->flag & BAM_FREAD1) ++(s)->n_read1[w];              \
			if ((c)->flag & BAM_FREAD2) ++(s)->n_read2[w];              \
			if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP)) ++(s)->n_sgltn[w];  \
			if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP)) { \
				++(s)->n_pair_map[w];                                   \
				if ((c)->mtid != (c)->tid) {                            \
					++(s)->n_diffchr[w];                                \
					if ((c)->qual >= 5) ++(s)->n_diffhigh[w];           \
				}                                                       \
			}                                                           \
		} else if (!((c)->flag & BAM_FUNMAP)) {                         \
			++(s)->n_u_mapped[w];                                       \
			if((c)->qual >= 30) ++(s)->n_u_map_good[w];                 \
        }                                                               \
		if (!((c)->flag & BAM_FUNMAP)) {                                \
			++(s)->n_mapped[w];                                         \
			if((c)->qual >= 30) ++(s)->n_map_good[w];                   \
		}                                                               \
		if ((c)->flag & BAM_FDUP) ++(s)->n_dup[w];                      \
	} while (0)

bam_flagstat_t *bam_flagstat_core(bamFile fp)
{
	bam_flagstat_t *s;
	bam1_t *b;
	bam1_core_t *c;
	int ret;
	s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));
	b = bam_init1();
	c = &b->core;
	while ((ret = bam_read1(fp, b)) >= 0)
		flagstat_loop(s, c);
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");
	return s;
}
int bam_flagstat(int argc, char *argv[], int vanilla)
{
	bamFile fp;
	bam_header_t *header;
	bam_flagstat_t *s;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s flagstat <in.bam>\n", invocation_name);
		return 1;
	}
	fp = strcmp(argv[1], "-")? bam_open(argv[1], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	s = bam_flagstat_core(fp);
	printf("%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
	printf("%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
	printf("%lld + %lld mapped (%.2f%%:%.2f%%)\n", s->n_mapped[0], s->n_mapped[1], (float)s->n_mapped[0] / s->n_reads[0] * 100.0, (float)s->n_mapped[1] / s->n_reads[1] * 100.0);
    if( !vanilla ) {
        printf("%lld + %lld mapped at MAPQ>=30 (%.2f%%:%.2f%%)\n", s->n_map_good[0], s->n_map_good[1], (float)s->n_map_good[0] / s->n_reads[0] * 100.0, (float)s->n_map_good[1] / s->n_reads[1] * 100.0);
        printf("%lld + %lld unpaired in sequencing\n", s->n_reads[0]-s->n_pair_all[0], s->n_reads[1]-s->n_pair_all[1]);
        printf("%lld + %lld unpaired and mapped (%.2f%%:%.2f%%)\n", s->n_u_mapped[0], s->n_u_mapped[1], (float)s->n_u_mapped[0] / (s->n_reads[0]-s->n_pair_all[0]) * 100.0, (float)s->n_u_mapped[1] / (s->n_reads[1]-s->n_pair_all[1]) * 100.0);
        printf("%lld + %lld unpaired and mapped at MAPQ>=30 (%.2f%%:%.2f%%)\n", s->n_u_map_good[0], s->n_u_map_good[1], (float)s->n_u_map_good[0] / (s->n_reads[0]-s->n_pair_all[0]) * 100.0, (float)s->n_u_map_good[1] / (s->n_reads[1]-s->n_pair_all[1]) * 100.0);
    }
	printf("%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
	printf("%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
	printf("%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
	printf("%lld + %lld properly paired (%.2f%%:%.2f%%)\n", s->n_pair_good[0], s->n_pair_good[1], (float)s->n_pair_good[0] / s->n_pair_all[0] * 100.0, (float)s->n_pair_good[1] / s->n_pair_all[1] * 100.0);
	printf("%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
	printf("%lld + %lld singletons (%.2f%%:%.2f%%)\n", s->n_sgltn[0], s->n_sgltn[1], (float)s->n_sgltn[0] / s->n_pair_all[0] * 100.0, (float)s->n_sgltn[1] / s->n_pair_all[1] * 100.0);
	printf("%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
	printf("%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
	free(s);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}

const char *labels[] = {
         "read group",
         "QC pass/fail",   
	     "total in this read group",
         "fraction in this read group",
	     "duplicates",
	     "mapped",
	     "fraction mapped",
	     "mapped at MAPQ>=30",
	     "fraction mapped at MAPQ>=30",
	     "unpaired in sequencing",
	     "unpaired in sequencing and mapped",
         "fraction of unpaired in sequencing mapped",
	     "unpaired in sequencing and mapped at MAPQ>=30",
	     "fraction of unpaired in sequencing mapped at MAPQ>=30",
	     "paired in sequencing",
	     "read1",
	     "read2",
	     "properly paired",
	     "fraction properly paired",
	     "paired with itself and mate mapped",
	     "singletons",
         "fraction singletons",
	     "paired with mate mapped to a different chr",
	     "paired with mate mapped to a different chr (mapQ>=5)" } ;

static const int nlabels = sizeof(labels)/sizeof(labels[0]); 

void print_header(FILE *hdl)
{
    int i;
    for(i=0;i!=nlabels;i++)
        fprintf(hdl, "# %2d - %s\n", i+1, labels[i]);
    fprintf(hdl, "# ");
    for(i=1;i!=nlabels;i++) fprintf(hdl, "%d\t",i);
    fprintf(hdl, "%d\n",i);
}

void print_flagstat( FILE* hdl, const char* rg, char* m, bam_flagstat_t *s, int w, int total )
{
    fprintf(hdl, "%s\t%s\t%lld\t%.2f%%\t%lld\t%lld\t%.2f%%\t%lld\t%.2f%%\t%lld\t%lld\t%.2f%%\t"
            "%lld\t%.2f%%\t%lld\t%lld\t%lld\t%lld\t%.2f%%\t%lld\t%lld\t%.2f%%\t%lld\t%lld\n",
            *rg ? rg : "\"\"", m, s->n_reads[w], (float)s->n_reads[w] / total * 100.0, s->n_dup[w], 
            s->n_mapped[w], (float)s->n_mapped[w] / s->n_reads[w] * 100.0,
            s->n_map_good[w], (float)s->n_map_good[w] / s->n_reads[w] * 100.0,
            s->n_reads[w]-s->n_pair_all[w],
            s->n_u_mapped[w], (float)s->n_u_mapped[w] / (s->n_reads[w]-s->n_pair_all[w]) * 100.0,
            s->n_u_map_good[w], (float)s->n_u_map_good[w] / (s->n_reads[w]-s->n_pair_all[w]) * 100.0,
            s->n_pair_all[w], s->n_read1[w], s->n_read2[w],
            s->n_pair_good[w], (float)s->n_pair_good[w] / s->n_pair_all[w] * 100.0,
            s->n_pair_map[w], s->n_sgltn[w], (float)s->n_sgltn[w] / s->n_pair_all[w] * 100.0,
            s->n_diffchr[w], s->n_diffhigh[w]) ;
}


void flagstatx_init( struct flagstatx_acc *f ) {
    memset( &f->total, 0, sizeof(bam_flagstat_t) );
    f->h = kh_init(rg2stat);
}

void flagstatx_destroy( struct flagstatx_acc *f ) {
    khint_t k;
	for (k = kh_begin(f->h); k != kh_end(f->h); ++k) 
        if(kh_exist(f->h,k)) 
            free(kh_key(f->h,k));
    kh_destroy(rg2stat, f->h);
}

void flagstatx_step( struct flagstatx_acc *f, const char* rg, bam1_t *b ) {
    int ret;
    khint_t k = kh_get(rg2stat, f->h, rg);
    if(k==kh_end(f->h)) { 
        k = kh_put(rg2stat,f->h,strdup(rg),&ret) ;
        memset(&kh_val(f->h,k), 0, sizeof(bam_flagstat_t));
    }
    bam_flagstat_t *s = &kh_val(f->h,k);
    flagstat_loop(s, &b->core);
    flagstat_loop(&f->total, &b->core);
}

void flagstatx_print( struct flagstatx_acc *f, FILE *hdl )
{
    khint_t k;
    int tot = f->total.n_reads[0] + f->total.n_reads[1];
    print_header(hdl);
	for (k = kh_begin(f->h); k != kh_end(f->h); ++k) {
        if(kh_exist(f->h,k)) {
            print_flagstat(hdl, kh_key(f->h,k), "pass", &kh_value(f->h,k), 0, tot);
            print_flagstat(hdl, kh_key(f->h,k), "fail", &kh_value(f->h,k), 1, tot);
        }
    }
    print_flagstat(hdl, "TOTAL", "pass", &f->total, 0, tot);
    print_flagstat(hdl, "TOTAL", "fail", &f->total, 1, tot);
}


int bam_flagstatx(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *header;
    struct flagstatx_acc flagstatx_acc ;
	bam1_t *b;
    int ret;

	if (argc == 1) {
		fprintf(stderr, "Usage: %s flagstatx <in.bam>\n", invocation_name);
		return 1;
	}

	fp = strcmp(argv[1], "-")? bam_open(argv[1], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);

    flagstatx_init( &flagstatx_acc ) ;

	b = bam_init1();
	while ((ret = bam_read1(fp, b)) >= 0) {
        flagstatx_step( &flagstatx_acc, get_rg(b), b ) ;
    }
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");

    flagstatx_print( &flagstatx_acc, stdout );
    flagstatx_destroy( &flagstatx_acc );
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}

static int usage(int lng) {
    fprintf(stderr, "Usage: %s covstat [-X|-L|-m N|-M N...] <in.bam>\n", invocation_name);
    fprintf(stderr, "  Counts coverage (with -X) or aligned length (with -L).\n");
    if(lng) {
        fprintf(stderr, "  Options are: -X        Count coverage in \"X\"\n");
        fprintf(stderr, "               -L        Count total aligned length\n");
        fprintf(stderr, "               -m N      Consider only reads of at least N bases\n");
        fprintf(stderr, "               -M N      Consider only reads of at most N bases\n");
    }
    return !lng;
}

void covstat_init( struct covstat_acc *c )
{
    c->h = kh_init(rg2cov);
    c->min_len=0;
    c->max_len=0x7fffffff;
    c->norm=1;
}

void covstat_step( struct covstat_acc *c, const char *rg, const bam_header_t *header, bam1_t *b ) {
    khint_t k;
    int ret;
    if (!(b->core.flag & BAM_FUNMAP) && b->core.tid >= 0 && b->core.tid < header->n_targets) {
        int len = b->core.flag & BAM_FPROPER_PAIR ? b->core.isize : b->core.l_qseq ;
        if( len >= c->min_len && len <= c->max_len ) {
            k = kh_get(rg2cov, c->h, rg);
            if(k==kh_end(c->h)) { 
                k = kh_put(rg2cov,c->h,strdup(rg),&ret) ;
                kh_val(c->h,k) = calloc( header->n_targets*2, sizeof(long long) );
            }
            kh_val(c->h,k)[ b->core.tid*2 + (b->core.flag & BAM_FQCFAIL ? 1 : 0) ] += 
                bam_calend(&b->core, bam1_cigar(b)) - b->core.pos ;
        }
    }
}

void covstat_print( struct covstat_acc *acc, FILE* f, bam_header_t *header ) {
    khint_t k;
    int i;
    fputs("# chr\tQC\t",f);
    for (k = kh_begin(acc->h); k != kh_end(acc->h); ++k) {
        if(kh_exist(acc->h,k)) {
            fputs(kh_key(acc->h,k),f);
            fputc('\t',f);
        }
    }
    fputs("TOTAL\n",f);

    for(i=0;i!=2*header->n_targets;++i) {
        long long t=0;
        fputs(header->target_name[i/2], f);
        fputs(i&1 ? "\tfail\t" : "\tpass\t", f);

        for (k = kh_begin(acc->h); k != kh_end(acc->h); ++k) {
            if(kh_exist(acc->h,k)) {
                if (acc->norm)
                    fprintf(f,"%.2f\t", (double)kh_value(acc->h,k)[i] / header->target_len[i/2]) ;
                else
                    fprintf(f,"%lld\t", kh_value(acc->h,k)[i]) ;
                t += kh_value(acc->h,k)[i] ;
            }
        }
        if (acc->norm)
            fprintf(f,"%.2f\n", (double)t / header->target_len[i/2]) ;
        else
            fprintf(f,"%lld\n", t) ;
    }
}

void covstat_destroy(struct covstat_acc *c) {
    khint_t k;
    for (k = kh_begin(c->h); k != kh_end(c->h); ++k)
        if(kh_exist(c->h,k)) {
            free(kh_key(c->h,k));
            free(kh_value(c->h,k));
        }
    kh_destroy(rg2cov, c->h);
}

int bam_covstat(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *header;
    int c;
    struct covstat_acc covstat_acc ;
    covstat_init( &covstat_acc ) ;

	while ((c = getopt(argc, argv, "m:M:LXh?")) >= 0) {
		switch (c) {
            case 'L': covstat_acc.norm=0; break;
            case 'X': covstat_acc.norm=1; break;
            case 'm': covstat_acc.min_len = atoi(optarg); break;
            case 'M': covstat_acc.max_len = atoi(optarg); break;
            case 'h':
            case '?': return usage(1);
            default: return usage(0);
        }
    }

	fp = optind != argc && strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);

	bam1_t *b;
	int ret;
	b = bam_init1();
	while ((ret = bam_read1(fp, b)) >= 0) {
        covstat_step( &covstat_acc, get_rg(b), header, b ) ;
    }
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");

    covstat_print(&covstat_acc, stdout, header);
    covstat_destroy(&covstat_acc);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}
