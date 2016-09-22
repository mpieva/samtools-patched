#include "bam.h"

struct index_acc {
	bam_index_t *idx;

	uint32_t last_bin, save_bin;
	int32_t last_coor, last_tid, save_tid;
	uint64_t save_off, last_off, n_mapped, n_unmapped, off_beg, off_end, n_no_coor;
};


void index_acc_init_A( struct index_acc *a ) ;
void index_acc_init_B( struct index_acc *a, int32_t n_targets, int64_t off0 ) ;
void index_acc_step_A( struct index_acc *acc, bam1_t *b ) ;
void index_acc_step( struct index_acc *acc, bam1_t *b, int64_t ) ;
bam_index_t *index_acc_finish( struct index_acc*, int64_t ) ;

void bam_index_save(const bam_index_t*, FILE*) ;
