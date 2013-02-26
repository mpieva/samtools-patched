#include <unistd.h>
#include <assert.h>
#include "bam.h"
#include "khash.h"

typedef struct {
	long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_u_mapped[2], n_u_map_good[2];
	long long n_sgltn[2], n_read1[2], n_read2[2];
	long long n_dup[2], n_map_good[2];
	long long n_diffchr[2], n_diffhigh[2];
} bam_flagstat_t;

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

KHASH_MAP_INIT_STR(rg2stat,  bam_flagstat_t)

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

void print_header()
{
    int i;
    for(i=0;i!=nlabels;i++)
        printf("# %2d - %s\n", i+1, labels[i]);
    printf("# ");
    for(i=1;i!=nlabels;i++) printf("%d\t",i);
    printf("%d\n",i);
}

void print_flagstat( const char* rg, char* m, bam_flagstat_t *s, int w, int total )
{
    printf("%s\t%s\t%lld\t%.2f%%\t%lld\t%lld\t%.2f%%\t%lld\t%.2f%%\t%lld\t%lld\t%.2f%%\t"
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


int bam_flagstatx(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *header;
    bam_flagstat_t total;
    memset( &total, 0, sizeof(bam_flagstat_t) );

    khash_t(rg2stat) *h = kh_init(rg2stat);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s flagstatx <in.bam>\n", invocation_name);
		return 1;
	}
	fp = strcmp(argv[1], "-")? bam_open(argv[1], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);

	bam1_t *b;
    khint_t k;
	int ret;
	b = bam_init1();
	while ((ret = bam_read1(fp, b)) >= 0) {
        char rg_[] = {0,0};
        char *rg = (char*)bam_aux_get(b,"RG");
        if(!rg) rg = rg_ ;
        else if(*rg=='A'||*rg=='c'||*rg=='C') {
            rg_[0]=rg[1];
            rg=rg_;
        }
        else if(*rg=='Z'||*rg=='H') ++rg;
        else rg = rg_ ;

        k = kh_get(rg2stat, h, rg);
        if(k==kh_end(h)) { 
            k = kh_put(rg2stat,h,strdup(rg),&ret) ;
            memset(&kh_val(h,k), 0, sizeof(bam_flagstat_t));
        }
        bam_flagstat_t *s = &kh_val(h,k);
		flagstat_loop(s, &b->core);
        flagstat_loop(&total, &b->core);
    }
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");

    int tot = total.n_reads[0]+total.n_reads[1];
    print_header();
	for (k = kh_begin(h); k != kh_end(h); ++k) {
        if(kh_exist(h,k)) {
            print_flagstat(kh_key(h,k), "pass", &kh_value(h,k), 0, tot);
            print_flagstat(kh_key(h,k), "fail", &kh_value(h,k), 1, tot);
        }
    }
    print_flagstat("TOTAL", "pass", &total, 0, tot);
    print_flagstat("TOTAL", "fail", &total, 1, tot);

	bam_header_destroy(header);
	bam_close(fp);
    kh_destroy(rg2stat, h);
	return 0;
}

/* Maps read group to raw coverage.  Each terget sequence gets two
 * entries in the coverage array, for quality pass and quality fail, in
 * that order. */
KHASH_MAP_INIT_STR(rg2cov, long long*);

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

int bam_covstat(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *header;
    khash_t(rg2cov) *h = kh_init(rg2cov);
    int c, norm=1, min_len=0, max_len=0x7fffffff;

	while ((c = getopt(argc, argv, "m:M:LXh?")) >= 0) {
		switch (c) {
            case 'L': norm=0; break;
            case 'X': norm=1; break;
            case 'm': min_len = atoi(optarg); break;
            case 'M': max_len = atoi(optarg); break;
            case 'h':
            case '?': return usage(1);
            default: return usage(0);
        }
    }

	fp = optind != argc && strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);

	bam1_t *b;
    khint_t k;
	int ret,i;
	b = bam_init1();
	while ((ret = bam_read1(fp, b)) >= 0) {
        if (!(b->core.flag & BAM_FUNMAP) && b->core.tid >= 0 && b->core.tid < header->n_targets) {
            int len = b->core.flag & BAM_FPROPER_PAIR ? b->core.isize : b->core.l_qseq ;
            if( len >= min_len && len <= max_len ) {
                char rg_[] = {0,0};
                char *rg = (char*)bam_aux_get(b,"RG");
                if(!rg) rg = rg_ ;
                else if(*rg=='A'||*rg=='c'||*rg=='C') {
                    rg_[0]=rg[1];
                    rg=rg_;
                }
                else if(*rg=='Z'||*rg=='H') ++rg;
                else rg = rg_ ;

                k = kh_get(rg2cov, h, rg);
                if(k==kh_end(h)) { 
                    k = kh_put(rg2cov,h,strdup(rg),&ret) ;
                    kh_val(h,k) = calloc( header->n_targets*2, sizeof(long long) );
                }
                kh_val(h,k)[ b->core.tid*2 + (b->core.flag & BAM_FQCFAIL ? 1 : 0) ] += 
                    bam_calend(&b->core, bam1_cigar(b)) - b->core.pos ;
            }
        }
    }
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");

    fputs("# chr\tQC\t",stdout);
    for (k = kh_begin(h); k != kh_end(h); ++k) {
        if(kh_exist(h,k)) {
            fputs(kh_key(h,k),stdout);
            fputc('\t',stdout);
        }
    }
    fputs("TOTAL\n",stdout);

    for(i=0;i!=2*header->n_targets;++i) {
        long long t=0;
        fputs(header->target_name[i/2], stdout);
        fputs(i&1 ? "\tfail\t" : "\tpass\t", stdout);

        for (k = kh_begin(h); k != kh_end(h); ++k) {
            if(kh_exist(h,k)) {
                if (norm)
                    printf("%.2f\t", (double)kh_value(h,k)[i] / header->target_len[i/2]) ;
                else
                    printf("%lld\t", kh_value(h,k)[i]) ;
                t += kh_value(h,k)[i] ;
            }
        }
        if (norm)
            printf("%.2f\n", (double)t / header->target_len[i/2]) ;
        else
            printf("%lld\n", t) ;
    }

	bam_header_destroy(header);
	bam_close(fp);

    for (k = kh_begin(h); k != kh_end(h); ++k)
        if(kh_exist(h,k))
            free(kh_value(h,k));
    kh_destroy(rg2cov, h);
	return 0;
}
