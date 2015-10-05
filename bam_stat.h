#include "bam.h"
#include "khash.h"

typedef struct {
	long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_u_mapped[2], n_u_map_good[2];
	long long n_sgltn[2], n_read1[2], n_read2[2];
	long long n_dup[2], n_map_good[2];
	long long n_diffchr[2], n_diffhigh[2];
} bam_flagstat_t;

KHASH_MAP_INIT_STR(rg2stat,  bam_flagstat_t)

struct flagstatx_acc {
    bam_flagstat_t total;
    khash_t(rg2stat) *h ;
} ;

void flagstatx_init( struct flagstatx_acc* );
void flagstatx_destroy( struct flagstatx_acc* );
void flagstatx_print( struct flagstatx_acc*, FILE* );
void flagstatx_step( struct flagstatx_acc*, const char* rg, bam1_t* );


/* Maps read group to raw coverage.  Each terget sequence gets two
 * entries in the coverage array, for quality pass and quality fail, in
 * that order. */
KHASH_MAP_INIT_STR(rg2cov, long long*);

struct covstat_acc {
    khash_t(rg2cov) *h;
    int min_len, max_len, n_targets, norm;
} ;

void covstat_init( struct covstat_acc* );
void covstat_step( struct covstat_acc*, const char *rg, const bam_header_t*, bam1_t* );
void covstat_print( struct covstat_acc*, FILE*, bam_header_t* );
void covstat_destroy( struct covstat_acc* );

inline static const char *get_rg( const bam1_t *b )
{
    static char rg[2];
    rg[0]=rg[1]=0;
    char *r = (char*)bam_aux_get(b,"RG");
    if(!r) return rg; 
    else if(*r=='A'||*r=='c'||*r=='C') {
        rg[0]=r[1];
        return rg;
    }
    else if(*r=='Z'||*r=='H') return r+1;
    else return rg;
}


