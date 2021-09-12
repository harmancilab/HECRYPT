#include <stdio.h>
#include <stdlib.h>
#include "hcrypt_mapped_read_tools.h"
#include <vector>
#include <ctype.h>
#include <math.h>
#include "hcrypt_signal_track_tools.h"
#include "hcrypt_annot_region_tools.h"
#include "hcrypt_genome_sequence_tools.h"
#include "hcrypt_variation_tools.h"
#include "hcrypt_genomics_coords.h"
#include "hcrypt_utils.h"
#include "hcrypt_nomenclature.h"
#include "hcrypt_rng.h"
#include "hcrypt_seed_manager.h"
#include "hcrypt_nucleotide.h"
#include <string.h>
#include <algorithm>
#include "hcrypt_ansi_string.h"

using namespace std; 

bool __DUMP_MAPPED_READ_TOOLS_MSGS__ = false;

FILE* get_processed_reads_ptr_wrapper(char* cur_chr_reads_fp)
{
	FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
//	if (t_string::compare_strings(cur_chr_reads_fp, "stdin"))
//	{
//		f_cur_chr_reads = stdin;
//	}
//	else if (t_string::ends_with(cur_chr_reads_fp, ".gz"))
//	{
//		char ungzip_cmd[1000];
//		sprintf(ungzip_cmd, "gzip -cd %s", cur_chr_reads_fp);
//#ifdef _WIN32
//		f_cur_chr_reads = _popen(ungzip_cmd, "r");
//#else 
//		f_cur_chr_reads = popen(ungzip_cmd, "r");
//#endif
//	}
//	else
//	{
//		f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
//	}

	return f_cur_chr_reads;
}


void count_10X_reads_per_cell_per_barcodes_list(char* TenX_SAM_file, char* per_cell_barcode_list_fp, char* op_fp)
{
	fprintf(stderr, "Loading and sorting barcodes from %s.\n", per_cell_barcode_list_fp);
	vector<char*>* per_cell_barcodes = buffer_file(per_cell_barcode_list_fp);
	sort(per_cell_barcodes->begin(), per_cell_barcodes->end(), t_string::sort_strings_per_prefix);
	int* n_reads_per_barcode = new int[per_cell_barcodes->size() + 2];
	memset(n_reads_per_barcode, 0, sizeof(int) * per_cell_barcodes->size());

	fprintf(stderr, "Generating the read count per cell for %d cell barcodes.\n", per_cell_barcodes->size());

	unsigned int n_reads_processed = 0;
	unsigned int n_unmatched_BC_reads = 0;
	FILE* f_sam = open_f(TenX_SAM_file, "r");
	while (1)
	{
		char* SAM_line = getline(f_sam);
		if (SAM_line == NULL)
		{
			break;
		}

		n_reads_processed++;
		if (n_reads_processed % (1000 * 1000) == 0)
		{
			fprintf(stderr, "Processing %d. read; %d unmatched BC reads.                       \r", n_reads_processed, n_unmatched_BC_reads);
		}

		// Parse the barcode and update counts.
		char barcode_buffer[100];
		if (get_10X_cell_barcode_per_SAM_read(SAM_line, barcode_buffer))
		{
			int cell_i = -1;
			int search_cell_i = t_string::fast_search_string_per_prefix(barcode_buffer, per_cell_barcodes, 0, per_cell_barcodes->size() - 1);
			while (search_cell_i > 0 &&
				(t_string::sort_strings_per_prefix(barcode_buffer, per_cell_barcodes->at(search_cell_i)) ||
					t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
			{
				search_cell_i--;
			} // search_cell_i loop.

			while (search_cell_i < per_cell_barcodes->size() &&
				(t_string::sort_strings_per_prefix(per_cell_barcodes->at(search_cell_i), barcode_buffer) ||
					t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
			{
				if (t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer))
				{
					cell_i = search_cell_i;
					break;
				}
				else
				{
					search_cell_i++;
				}
			} // search_cell_i loop.

			  // If we found a cell, update the count of it.
			if (cell_i != -1)
			{
				n_reads_per_barcode[cell_i]++;
			}
			else
			{
				n_unmatched_BC_reads++;
			}
		} // tag parsing check.

		delete[] SAM_line;
	} // file reading loop.
	close_f(f_sam, TenX_SAM_file);

	// Write output.
	fprintf(stderr, "Saving the read counts per barcode to %s\n", op_fp);
	FILE* f_per_BC_read_counts = open_f(op_fp, "w");
	for (int bc_i = 0; bc_i < per_cell_barcodes->size(); bc_i++)
	{
		fprintf(f_per_BC_read_counts, "%s\t%d\n", per_cell_barcodes->at(bc_i), n_reads_per_barcode[bc_i]);
	} // bc_i loop.
	fclose(f_per_BC_read_counts);
}

bool get_10X_cell_barcode_per_SAM_read(char* TenX_SAM_read_line, char* barcode_buffer)
{
	char temp_tag_entry_buffer[1000];
	bool parsed_BC_tag = get_SAM_read_tag_entry(TenX_SAM_read_line, temp_tag_entry_buffer, "CB:Z:");

	// Copy the barcode.
	if (parsed_BC_tag)
	{
		strcpy(barcode_buffer, &(temp_tag_entry_buffer[5]));
	}

	return(parsed_BC_tag);
}


bool get_SAM_read_tag_entry(char* sam_read_line, char* tag_entry_buffer, const char* tag_prefix)
{
	int i_cur_char = 0;
	char cur_SAM_entry[1010];
	int l_entry_buff = 1000;

	// Find 11th token.
	int i_tok = 0;
	while (i_tok < 11)
	{
		t_string::get_next_token(sam_read_line, NULL, l_entry_buff, "\t", i_cur_char);
		i_tok++;
	} // i_tok loop.

	// Loop through next set of tokens.
	bool found_tag = false;
	while (t_string::get_next_token(sam_read_line, cur_SAM_entry, l_entry_buff, "\t", i_cur_char))
	{
		// Copy the barcode sequence.
		if (t_string::starts_with(cur_SAM_entry, tag_prefix))
		{
			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "%s tag: %s\n", tag_prefix, cur_SAM_entry);
			}

			strcpy(tag_entry_buffer, cur_SAM_entry);
			found_tag = true;
			break;
		}
	} // i_tok loop. 

	return(found_tag);
}

void close_processed_reads_ptr_wrapper(FILE* f_cur_chr_reads, char* cur_chr_reads_fp)
{
	if (t_string::compare_strings(cur_chr_reads_fp, "stdin"))
	{
	}
	else if (t_string::ends_with(cur_chr_reads_fp, ".gz"))
	{
#ifdef _WIN32
		_pclose(f_cur_chr_reads);
#else 
		pclose(f_cur_chr_reads);
#endif
	}
	else
	{
		fclose(f_cur_chr_reads);
	}
}

bool set_preprocessed_read_file_path_per_dir_chr(char* preprocessed_reads_dir, char* chrom, char* preprocessed_reads_fp)
{
	sprintf(preprocessed_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chrom);

	if (check_file(preprocessed_reads_fp))
	{
		return true;
	}
	
	sprintf(preprocessed_reads_fp, "%s/%s_mapped_reads.txt.gz", preprocessed_reads_dir, chrom);

	if (check_file(preprocessed_reads_fp))
	{
		return true;
	}

	return false;
}

vector<int>* count_preprocessed_reads(char* preprocessed_reads_dir, vector<char*>* _chr_ids)
{
	// Load all the reads per chromosome: Load all the chromosomes.
	char chr_ids_list_fp[1000];
	sprintf(chr_ids_list_fp, "%s/chr_ids.txt", preprocessed_reads_dir);

	// Load the chromosome id's.
	vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
	if (chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_list_fp);
		return(NULL);
	}

	vector<int>* n_reads_per_chromosome = new vector<int>();

	t_rng* sampling_rng = new t_rng(t_seed_manager::seed_me());

	int n_total_reads = 0;
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		_chr_ids->push_back(t_string::copy_me_str(chr_ids->at(i_chr)));

		char cur_chr_reads_fp[1000];
		//sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(preprocessed_reads_dir, chr_ids->at(i_chr), cur_chr_reads_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", preprocessed_reads_dir, chr_ids->at(i_chr));
			exit(0);
		}

		fprintf(stderr, "Counting the reads in %s (%s)\n", chr_ids->at(i_chr), cur_chr_reads_fp);

		//int n_cur_chr_reads = 0;
		FILE* f_cur_chr_reads = get_processed_reads_ptr_wrapper(cur_chr_reads_fp);

		int n_cur_chr_reads = get_n_lines(f_cur_chr_reads);

		close_processed_reads_ptr_wrapper(f_cur_chr_reads, cur_chr_reads_fp);

		n_reads_per_chromosome->push_back(n_cur_chr_reads);

		// Update total # of reads.
		n_total_reads += n_cur_chr_reads;

		//fprintf(stderr, "%s: %d reads. (%d reads in total)\n", chr_ids->at(i_chr), n_cur_chr_reads, n_total_reads);
	} // i_chr loop.

	return(n_reads_per_chromosome);
}

void get_signal_statistics_in_regions_per_pileups_dir(char* chrom_info_fp, char* pileups_dir_path, char* regs_bed_fp)
{
	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chrom_info_fp, chr_ids, chr_lengths);
	fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

	vector<t_annot_region*>* regs = load_BED(regs_bed_fp);
	fprintf(stderr, "Loaded %d regions.\n", regs->size());

	double total_signal = 0;
	double total_signal_in_regions = 0;
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_regs = get_regions_per_chromosome(regs, chr_ids->at(i_chr));

		char pileup_fp[1000];
		sprintf(pileup_fp, "%s/%s_allele_counts.bin", pileups_dir_path, chr_ids->at(i_chr));
		int* covg = NULL;
		int l_sig = 0;
		if (!check_file(pileup_fp))
		{
			sprintf(pileup_fp, "%s/%s_allele_counts.bin.gz", pileups_dir_path, chr_ids->at(i_chr));

			if (!check_file(pileup_fp))
			{
				fprintf(stderr, "Could not find the pileup @ %s\n", pileup_fp);
				exit(0);
			}
		}

		covg = load_coverage_per_compressed_pileup_file(pileup_fp, l_sig);
		fprintf(stderr, "Loaded coverage of %d nucleotides.\n", l_sig);

		for (int i = 0; i <= l_sig; i++)
		{
			total_signal += covg[i];
		} // i loop.

		for (int i_reg = 0; i_reg < cur_chr_regs->size(); i_reg++)
		{
			for (int i = cur_chr_regs->at(i_reg)->start; i <= cur_chr_regs->at(i_reg)->end; i++)
			{
				total_signal_in_regions += covg[i];
			} // i loop.
		} // i_reg loop.
	} // i_chr loop.

	fprintf(stderr, "%d regions, %lf / %lf = %lf signal in regions.\n", regs->size(), total_signal_in_regions, total_signal, total_signal_in_regions / total_signal);
}

bool sort_read_line_entries_per_id(t_read_line_w_id* read1, t_read_line_w_id* read2)
{
	return(t_string::sort_strings(read1->id, read2->id));
}

#define __UCHAR_MAPPABILITY__

double* load_normalized_multimappability_profile(char* mapability_signal_profile_fp, int& l_mapability_profile)
{
	double* mapability_signal_profile = NULL;

	if(!check_file(mapability_signal_profile_fp))
	{
		l_mapability_profile = 0;
		return(NULL);
	}

#ifdef __DOUBLE_MAPPABILITY__
	// Load the mapability map signal profile, do filtering.
	mapability_signal_profile = load_per_nucleotide_binary_profile(mapability_signal_profile_fp, l_mapability_profile);

	// Do mapability aware median filtering on the current signal profile.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
	fprintf(stderr, "Scaling the mapability map with %d.\n", l_read_mapability_signal * 2);

	int mapability_scaling = l_read * 2;
	for(int i = 1; i <= l_mapability_profile; i++)
	{
		mapability_signal_profile[i] /= mapability_scaling;
	} // i loop.
#elif defined(__UCHAR_MAPPABILITY__)
	// Following loads the mappability signal profile from the char version of the multi-mappability profile.
	// Load the mapability map signal profile, do filtering.
	unsigned char* mapability_signal_char_profile = load_per_nucleotide_binary_uchar_profile(mapability_signal_profile_fp, l_mapability_profile);
	mapability_signal_profile = new double[l_mapability_profile + 2];
	for(int i = 1; i <= l_mapability_profile; i++)
	{
		unsigned char unsigned_char_val = (unsigned char)(mapability_signal_char_profile[i]);
		mapability_signal_profile[i] = (double)(unsigned_char_val);
		mapability_signal_profile[i] /= 100;

		if(mapability_signal_profile[i] < 0)
		{
			fprintf(stderr, "Sanity check failed.\n");
			exit(0);
		}
	} // i loop.
	delete [] mapability_signal_char_profile;
#else
	#error "Must define the type of mappability."
#endif

	return(mapability_signal_profile);
}

bool sort_read_lines(char* read1, char* read2)
{
	return(t_string::sort_strings(read1, read2));
}

#define MAX(x, y) (((x)>(y))?(x):(y))
#define MIN(x, y) (((x)<(y))?(x):(y))

int get_chr_i_per_group_entry(int self_flag_w_chr_w_mate_chr_val, bool self)
{
	if(self)
	{
		int flagged_val = (self_flag_w_chr_w_mate_chr_val & 0xFF00);
		int self_chr_i = flagged_val >> 8;
		return(self_chr_i);
	}
	else
	{
		int flagged_val = (self_flag_w_chr_w_mate_chr_val & 0xFF);
		int mate_chr_i = flagged_val;
		return(mate_chr_i);
	}
}

unsigned int get_flag_per_group_entry(int self_flag_w_chr_w_mate_chr_val)
{
	int flagged_val = (self_flag_w_chr_w_mate_chr_val & 0xFFFF0000);
	unsigned int flag_val = flagged_val >> 16;
	return(flag_val);
}

bool sort_barcode_entries_per_barcode(t_10X_BX_barcode_group_entry* ent1, t_10X_BX_barcode_group_entry* ent2)
{
	return(ent1->BX_barcode < ent2->BX_barcode);
}

int recursive_locate_barcode_entry_index(double cur_barcode, 
										vector<t_10X_BX_barcode_group_entry*>* barcode_entries, 
										int start_i, int end_i)
{
	int mid_i = (start_i + end_i) / 2;

	if(mid_i == start_i || mid_i == end_i)
	{
		return(mid_i);
	}

	if(cur_barcode > barcode_entries->at(mid_i)->BX_barcode)
	{
		return recursive_locate_barcode_entry_index(cur_barcode, 
													barcode_entries, 
													mid_i, end_i);
	}
	else if(cur_barcode < barcode_entries->at(mid_i)->BX_barcode)
	{
		return recursive_locate_barcode_entry_index(cur_barcode, 
													barcode_entries, 
													start_i, mid_i);
	}
	else
	{
		return mid_i;
	}

}

double get_BX_barcode_val_per_barcode_str(char* barcode_str, int* per_char_val, double base_val, int& n_processed_chars)
{
	double cur_barcode = 0;
	double cur_base = 1;
	int i = 0;
	while(barcode_str[i] &&
		per_char_val[(int)(barcode_str[i])] < 4)
	{
		int num_val = per_char_val[(int)(barcode_str[i])];
		cur_barcode += (cur_base * num_val);

		cur_base *= base_val;

		i++;
	} // i loop.

	n_processed_chars = i;

	return(cur_barcode);
}


void get_barcode_string_per_value(double barcode_val, double base_val, char* barcode_str)
{
	char nuc_per_val[] = "ACGTN";

	double cur_val = barcode_val;
	int nuc_i = 0;
	while(cur_val > 0)
	{
		int cur_nuc_val = cur_val - floor(cur_val / base_val) * base_val;
		barcode_str[nuc_i] = nuc_per_val[cur_nuc_val];
		cur_val = floor(cur_val / base_val);

		nuc_i++;
	} // cur_val loop.
	barcode_str[nuc_i] = 0;

	//fprintf(stderr, "Uninverted: %s\n", barcode_str);
	//t_string::revert(barcode_str);
}

void get_unique_10X_BX_tags(char* tenX_sam_fp, int min_mapQ, char* op_fp)
{
	fprintf(stderr, "Parsing BX tags from 10X reads in %s using minimum mapping quality of %d\n", tenX_sam_fp, min_mapQ);

	FILE* f_tenX_sam = NULL;

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
		f_tenX_sam = stdin;
	}
	else
	{
		f_tenX_sam = open_f(tenX_sam_fp, "r");
	}

	int* per_char_val = get_per_char_as_nuc_2_num_coding_array();

	int n_reads_no_BX = 0;
	int n_processed_reads = 0;

	int l_barcode = 16;

	// This is the list of barcode groups.
	vector<double>* bx_tag_vals = new vector<double>();
	while(1)
	{
		char* cur_line = getline(f_tenX_sam);
		if(cur_line == NULL)
		{
			break;
		}

		n_processed_reads++;

		if(n_processed_reads % (1000*1000) == 0)
		{
			fprintf(stderr, "Processing %d. read. (%d w BX, %d w/o BX)                    \r", n_processed_reads, bx_tag_vals->size(), n_reads_no_BX);
		}

		int flag = 0;
		char read_id[1000];
		char flag_str[100];
		int _chr_index = 0;
		char _chr_index_str[100];
		int _mate_chr_index = 0;
		char _mate_chr_index_str[100];
		char _mapQ_str[100];
		char fragment[100000];
		char phred_quality_str[100000];
		char cigar_str[1000];
		char chrom[100];
		char mate_chrom[100];

		// Parse the SAM line.
		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %[^\t]", 
							read_id, flag_str, chrom, _chr_index_str, _mapQ_str, cigar_str, 
							mate_chrom, _mate_chr_index_str, 
							fragment, phred_quality_str) == 10)
		{
			int mapQ = atoi(_mapQ_str);

			if(mapQ < min_mapQ)
			{

			}
			else
			{
				// Parse this read.
				bool found_BX = false;

				char cur_entry[1000];
				int i_cur_char = 0;
				while(t_string::get_next_token(cur_line, cur_entry, 1000, "\t", i_cur_char))
				{
					if(cur_entry[0] == 'B' && t_string::starts_with(cur_entry, "BX:Z:"))
					{			
						// Create and search the current list of GEM groups.
						char* cur_full_barcode_str = cur_entry;
						int n_proc_chars = 0;
						double cur_barcode = get_BX_barcode_val_per_barcode_str(&(cur_full_barcode_str[5]), per_char_val, 4, n_proc_chars);
						if(n_proc_chars != l_barcode)
						{
							fprintf(stderr, "Could not read %d characters from barcode: %s\n", l_barcode, &(cur_full_barcode_str[5]));
						}

						//char decoded_str[1000];
						//get_barcode_string_per_value(cur_barcode, 4, decoded_str);

						//if(!t_string::compare_strings(decoded_str, &(cur_full_barcode_str[5])))
						//{
						//	fprintf(stderr, "Decoded does not match original:\n%s\n%s\n", 
						//			decoded_str, &(cur_full_barcode_str[5]));
						//	exit(0);
						//}

						// Find this barcode among existing barcodes.
						bx_tag_vals->push_back(cur_barcode);
						found_BX = true;
						break;
					}
				} // token loop.

				if(!found_BX)
				{
					//fprintf(stderr, "Could not parse GEM group entry: %s\n", cur_line);
					n_reads_no_BX++;
				}
			} // mapQ check.
		} // Read info parse check.

		// Free memory.
		delete [] cur_line;
	} // file reading loop.

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_tenX_sam);
	}

	fprintf(stderr, "%d reads with BX tags (%d reads missing BX), saving to %s\n", bx_tag_vals->size(), n_reads_no_BX, op_fp);

	// Dump the unique tag values.
	sort(bx_tag_vals->begin(), bx_tag_vals->end());

	FILE* f_op = open_f(op_fp, "w");
	double cur_tag_val = -1;
	for(int i_tag = 0; i_tag < bx_tag_vals->size(); i_tag++)
	{
		if(cur_tag_val != bx_tag_vals->at(i_tag))
		{
			fprintf(f_op, "%.1f\n", bx_tag_vals->at(i_tag));			
		}

		cur_tag_val = bx_tag_vals->at(i_tag);
	} // i_tag loop.
	fclose(f_op);
}

// Parse the BX tag::https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam
void parse_10x_linked_reads_per_BX_tag(char* tenX_sam_fp, vector<char*>* chr_ids, char* GEM_barcodes_list_fp, int min_mapQ, char* op_fp)
{
	fprintf(stderr, "Parsing linked reads per BX tag from 10X sam file @ %s (min mapQ %d) using %d chromosomes with barcodes in file %s.\n", tenX_sam_fp, min_mapQ, chr_ids->size(), GEM_barcodes_list_fp);

	FILE* f_tenX_sam = NULL;

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
		f_tenX_sam = stdin;
	}
	else
	{
		f_tenX_sam = open_f(tenX_sam_fp, "r");
	}

	int* per_char_val = get_per_char_as_nuc_2_num_coding_array();

	// Load the list of barcode groups.
	vector<t_10X_BX_barcode_group_entry*>* barcode_entries = new vector<t_10X_BX_barcode_group_entry*>();
	fprintf(stderr, "Loading GEM/BX barcodes.\n");
	vector<char*>* barcodes = buffer_file(GEM_barcodes_list_fp);
	if(barcodes == NULL)
	{
		fprintf(stderr, "Could not load bardcodes from %s\n", GEM_barcodes_list_fp);
		exit(0);
	}

	for(int ibc = 0; ibc < barcodes->size(); ibc++)
	{
		t_10X_BX_barcode_group_entry* new_entry = new t_10X_BX_barcode_group_entry();
		new_entry->BX_barcode = atof(barcodes->at(ibc));
		new_entry->self_chr_w_mate_chr_w_flag = new vector<int>();
		new_entry->posn = new vector<int>();
		new_entry->mate_posn = new vector<int>();

		barcode_entries->push_back(new_entry);
	} // ibc loop.
	fprintf(stderr, "Loaded %d entries.\n", barcode_entries->size());

	// Sort the entries.
	sort(barcode_entries->begin(), barcode_entries->end(), sort_barcode_entries_per_barcode);

	int l_barcode = 16;

	int n_processed_reads = 0;
	fprintf(stderr, "Loading and parsing 10X SAM file.\n");
	while(1)
	{
		char* cur_line = getline(f_tenX_sam);
		if(cur_line == NULL)
		{
			break;
		}

		n_processed_reads++;

		if(n_processed_reads % 1000000 == 0)
		{
			fprintf(stderr, "Processing %d. read: %d barcode entries                  \r", n_processed_reads, barcode_entries->size());
		}

		// Parse this read.
		int flag = 0;
		char read_id[1000];
		char flag_str[100];
		int _chr_index = 0;
		char _chr_index_str[100];
		int _mate_chr_index = 0;
		char _mate_chr_index_str[100];
		char _mapQ_str[100];
		char fragment[100000];
		char phred_quality_str[100000];
		char cigar_str[1000];
		char chrom[100];
		char mate_chrom[100];

		unsigned char cur_chrom_i = chr_ids->size();
		unsigned char mate_chrom_i = chr_ids->size();		

		// Parse the SAM line.
		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %[^\t]", 
							read_id, flag_str, chrom, _chr_index_str, _mapQ_str, cigar_str, 
							mate_chrom, _mate_chr_index_str, 
							fragment, phred_quality_str) == 10)
		{
			int mapQ = atoi(_mapQ_str);

			if(mapQ >= min_mapQ)
			{
				_chr_index = atoi(_chr_index_str);
				_mate_chr_index = atoi(_mate_chr_index_str);
				flag = atoi(flag_str);		

				// Translate the SAM indexing.
				_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);
				_mate_chr_index = translate_coord(_mate_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

				// Check the flag and determine the strand.
				char strand_char = 'F';
				if(flag & 0x10)
				{
					strand_char = 'R';
				}

				normalize_chr_id(chrom);							

				if(mate_chrom[0] == '=')
				{
					strcpy(mate_chrom, chrom);
				}
				else
				{
					normalize_chr_id(mate_chrom);
				}

				unsigned char mate_chrom_i = chr_ids->size();
				mate_chrom_i = t_string::get_i_str(chr_ids, mate_chrom);

				// Sanity check. Is this fragment mapped?
				unsigned char cur_chrom_i = t_string::get_i_str(chr_ids, chrom);

				// Make sure this read passes the quality checks.
				if((flag & 0x400) == 0 &&					// Optical OR PCR duplicate.
					(flag & 0x200) == 0 &&					// Vendor quality checks.
					(flag & 0x4) == 0 &&					// Read mapped
					cur_chrom_i != chr_ids->size())
				{
					// Find the group that this read belongs to.
					bool found_BX = false;

					// Find the "BX:Z:" substring index.
					int char_i = 0;
					while(cur_line[char_i])
					{
						if(cur_line[char_i] == 'B' &&
							cur_line[char_i+1] == 'X' &&
							cur_line[char_i+2] == ':' &&
							cur_line[char_i+3] == 'Z' &&
							cur_line[char_i+4] == ':')
						{
							found_BX = true;

							//// Find the next tab and set it to 0.
							//for(int char_i2 = char_i+5; 
							//	cur_line[char_i2] != 0; 
							//	char_i2++)
							//{
							//	if(cur_line[char_i2] == '\t')
							//	{
							//		cur_line[char_i2] = 0;
							//	}
							//} // char_i2 loop.
							
							break;
						} // BX:Z: check.
						char_i++;
					} // char_i loop.

					if(found_BX)
					{
						// Create and search the current list of GEM groups.
						// Following fixes the naming for the barcode.
						char* cur_full_barcode_str = &(cur_line[char_i]);
						char* cur_barcode_str = &(cur_full_barcode_str[5]);

						int n_proc_chars = 0;
						double cur_barcode = get_BX_barcode_val_per_barcode_str(cur_barcode_str, per_char_val, 4, n_proc_chars);
						if(n_proc_chars != l_barcode)
						{
							fprintf(stderr, "Could not read %d (%d) characters from barcode: %s\n", l_barcode, n_proc_chars, cur_barcode_str);
							exit(0);
						}

						// Find this barcode among existing barcodes.
						int ent_i = recursive_locate_barcode_entry_index(cur_barcode, barcode_entries, 0, barcode_entries->size() - 1);

						while(ent_i > 0 && 
							ent_i < barcode_entries->size() && 
							barcode_entries->at(ent_i)->BX_barcode >= cur_barcode)
						{
							ent_i--;
						} // ent_i rewind loop.

						bool found_matching_barcode_entry = false;
						while(ent_i >= 0 && 
							ent_i < barcode_entries->size() && 
							barcode_entries->at(ent_i)->BX_barcode <= cur_barcode)
						{
							if(barcode_entries->at(ent_i)->BX_barcode == cur_barcode)
							{
								found_matching_barcode_entry = true;

								// Pack some values.
								int self_flag_w_chr_w_mate_chr_val = 0;
								self_flag_w_chr_w_mate_chr_val = (flag << 16) | (cur_chrom_i << 8) | (mate_chrom_i);

								barcode_entries->at(ent_i)->self_chr_w_mate_chr_w_flag->push_back(self_flag_w_chr_w_mate_chr_val);
								barcode_entries->at(ent_i)->posn->push_back(_chr_index);
								barcode_entries->at(ent_i)->mate_posn->push_back(_mate_chr_index);

								// Test if we can successfully retrieve the values.
								int rec_flag = get_flag_per_group_entry(self_flag_w_chr_w_mate_chr_val);
								int rec_self_chr_i = get_chr_i_per_group_entry(self_flag_w_chr_w_mate_chr_val, true);
								int rec_mate_chr_i = get_chr_i_per_group_entry(self_flag_w_chr_w_mate_chr_val, false);
								if(rec_flag != flag)
								{
									fprintf(stderr, "Recovered flag does not match: %d, %d: %s\n", rec_flag, (int)flag, cur_line);
									exit(0);
								}

								if(rec_self_chr_i != cur_chrom_i ||
									rec_mate_chr_i != mate_chrom_i)
								{
									fprintf(stderr, "Recovered chromosome indices do not match: %d, %d; %d, %d\n%s\n", 
										(int)rec_self_chr_i, (int)cur_chrom_i, 
										(int)rec_mate_chr_i, (int)mate_chrom_i,
										cur_line);
									exit(0);
								}

								// Break the barcode search loop.
								break;
							} // barcode match check.

							ent_i++;
						} // ent_i rewind loop.

						if(found_matching_barcode_entry == false)
						{
							fprintf(stderr, "Could not find entry with barcode %lf\n", cur_barcode);
						}
					} // BX:Z check.

					if(!found_BX)
					{
						//fprintf(stderr, "Could not parse GEM group entry: %s\n", cur_line);
					}
				} // read unmapping check from flag.
			} // Minimum mapQ check.
		} // sam file parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	// Save the linked read information.	

	// Dump the text output.
#define __DUMP_TEXT__
#undef __DUMP_BINARY__

	FILE* f_op = NULL;
#ifdef __DUMP_TEXT__
	fprintf(stderr, "Saving text linked read information to %s\n", op_fp);
	f_op = open_f(op_fp, "w");
#elif defined(__DUMP_BINARY__)
	fprintf(stderr, "Saving binary linked read information to %s\n", op_fp);
	f_op = open_f(op_fp, "wb");
#endif 

	for(int bc_i = 0; bc_i < barcode_entries->size(); bc_i++)
	{
		if(barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size() <= 1)
		{
			continue;
		}

#ifdef __DUMP_TEXT__
		fprintf(f_op, "%.1f\t%d", barcode_entries->at(bc_i)->BX_barcode, (int)barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size());
#elif defined(__DUMP_BINARY__)
		// Write the barcode.
		fwrite(&(barcode_entries->at(bc_i)->BX_barcode), sizeof(double), 1, f_op);

		// Write the number of entries in the barcode group.
		int n_reads = (int)(barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size());
		fwrite(&(n_reads), sizeof(int), 1, f_op);
#endif 
		for(int i_r = 0; 
			i_r < (int)barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->size();
			i_r++)
		{
#ifdef __DUMP_TEXT__
			fprintf(f_op, "\t%d\t%d\t%d", 
					barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->at(i_r), 
					barcode_entries->at(bc_i)->posn->at(i_r), 
					barcode_entries->at(bc_i)->mate_posn->at(i_r));
#elif defined(__DUMP_BINARY__)
			int cur_self_chr_w_mate_chr_w_flag = barcode_entries->at(bc_i)->self_chr_w_mate_chr_w_flag->at(i_r);
			int cur_posn = barcode_entries->at(bc_i)->posn->at(i_r);
			int cur_mate_posn = barcode_entries->at(bc_i)->mate_posn->at(i_r);

			// Write the information about this entry.
			fwrite(&(cur_self_chr_w_mate_chr_w_flag), sizeof(int), 1, f_op);
			fwrite(&(cur_posn), sizeof(int), 1, f_op);
			fwrite(&(cur_mate_posn), sizeof(int), 1, f_op);
#endif 
		} // i_r loop.

#ifdef __DUMP_TEXT__
		fprintf(f_op, "\n");
#endif
	} // bc_i loop.
	fclose(f_op);

	if(t_string::compare_strings(tenX_sam_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_tenX_sam);
	}
} // parse_10x_linked_reads_per_BX_tag

vector<t_annot_region*>* load_per_chrom_per_bin_discordancy_stats(char* per_bin_discordancy_fp, 
																					vector<char*>* chr_ids,
																					int l_bin)
{
	vector<t_annot_region*>* per_bin_discordancy_stat_regs = new vector<t_annot_region*>();

	FILE* f_per_bin_discordancy = open_f(per_bin_discordancy_fp, "r");
	while(1)
	{
		char* cur_line = getline(f_per_bin_discordancy);
		if(cur_line == NULL)
		{
			break;
		}

		char cur_chr_id[1000];
		int cur_bin_start = 0;
		double avg_multimapp = 0;
		int n_cis_discordant = 0;
		int n_trans_discordant = 0;
		if(sscanf(cur_line, "%s %d %lf %d %d", cur_chr_id, &cur_bin_start, &avg_multimapp, &n_cis_discordant, &n_trans_discordant) != 5)
		{
			fprintf(stderr, "Could not parse.\n");
			exit(0);
		}

		t_per_bin_10x_discordancy_stat* discordant_stats = new t_per_bin_10x_discordancy_stat();
		discordant_stats->cur_bin_start = cur_bin_start;		
		discordant_stats->avg_multimapp = avg_multimapp;
		discordant_stats->n_cis_discordant_reads = n_cis_discordant;
		discordant_stats->n_trans_discordant_reads = n_trans_discordant;

		t_annot_region* cur_reg = get_empty_region();
		cur_reg->chrom = t_string::copy_me_str(cur_chr_id);

		// We put a padding to get around intersecting of neighboring regions.
		cur_reg->start = cur_bin_start + 10;
		cur_reg->end = cur_bin_start + l_bin - 1 - 10;
		cur_reg->strand = '+';
		cur_reg->data = discordant_stats;
		per_bin_discordancy_stat_regs->push_back(cur_reg);
	} // file reading loop.
	fclose(f_per_bin_discordancy);

	return(per_bin_discordancy_stat_regs);
}

void get_differential_discordancy_stats_per_10x_linked_reads(vector<char*>* chr_ids,																
																char* tumor_discordancy_fp, 
																char* normal_discordancy_fp,
																int l_bin,
																char* op_fp)
{
	vector<t_annot_region*>* tumor_discordant_stat_regs = load_per_chrom_per_bin_discordancy_stats(tumor_discordancy_fp, chr_ids, l_bin);
	fprintf(stderr, "Loaded %d tumor discordancy regions.\n", tumor_discordant_stat_regs->size());

	vector<t_annot_region*>* normal_discordant_stat_regs = load_per_chrom_per_bin_discordancy_stats(normal_discordancy_fp, chr_ids, l_bin);
	fprintf(stderr, "Loaded %d normal discordancy regions.\n", normal_discordant_stat_regs->size());

	fprintf(stderr, "Intersecting discordancy stat regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(tumor_discordant_stat_regs, normal_discordant_stat_regs, true);

	fprintf(stderr, "Processing %d intersects.\n", intersects->size());

	FILE* f_op = open_f(op_fp, "w");

	fprintf(f_op, "#CHROMOSOME\tBIN_START\t\
AVG_MULTIMAP\t\
#_TUMOR_CIS_DISCORDANT\t#_TUMOR_TRANS_DISCORDANT\t\
#_NORMAL_CIS_DISCORDANT\t#_NORMAL_TRANS_DISCORDANT\n");

	for(int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_annot_region* tumor_reg = ((t_intersect_info*)(intersects->at(i_int)->data))->src_reg;
		t_annot_region* normal_reg = ((t_intersect_info*)(intersects->at(i_int)->data))->dest_reg;
		int l_overlap = ((t_intersect_info*)(intersects->at(i_int)->data))->l_overlap;

		if(l_overlap > (l_bin / 2))
		{
			t_per_bin_10x_discordancy_stat* tumor_stat = (t_per_bin_10x_discordancy_stat*)(tumor_reg->data);
			t_per_bin_10x_discordancy_stat* normal_stat = (t_per_bin_10x_discordancy_stat*)(normal_reg->data);

			// Write the overlap; if there are more discordant reads in tumor than normal.
			if(tumor_stat->n_cis_discordant_reads > normal_stat->n_cis_discordant_reads ||
				tumor_stat->n_trans_discordant_reads > normal_stat->n_trans_discordant_reads)
			{
				fprintf(f_op, "%s\t%d\t%.3f\t%d\t%d\t%d\t%d\n", tumor_reg->chrom, tumor_reg->start, 
						tumor_stat->avg_multimapp,
						tumor_stat->n_cis_discordant_reads, tumor_stat->n_trans_discordant_reads,
						normal_stat->n_cis_discordant_reads, normal_stat->n_trans_discordant_reads);
			}
		} // overlap check.
	} // i_int loop.
	fclose(f_op);
}

bool sort_10x_linked_read_info_per_position(t_10x_linked_read_info* read1, t_10x_linked_read_info* read2)
{
	return(read1->posn < read2->posn);
}

void compute_discordancy_stats_per_10x_linked_reads(char* parsed_linked_reads_fp, 
													char* chr_ids_lengths_list_fp,
													int l_bin,
													char* multimapp_dir,
													char* op_fp)
{
	// Parameters for defining discordancy over linked read clusters.
	int min_supporting_l_broken_fragment = 250;
	int min_n_supporting_reads_in_broken_fragment = 3;
	int min_l_discordant_interread_dist_on_molecule = 1000;
	double min_mult_d12_per_distance = 1.5;

	// Load the id's and lengths.
	if(!check_file(chr_ids_lengths_list_fp))
	{
		fprintf(stderr, "Could not load chromosome ids/lengths from %s\n", chr_ids_lengths_list_fp);
		exit(0);
	}

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_list_fp, chr_ids, chr_lengths);

	// Generate the genomewide bins.
	vector<t_per_bin_10x_discordancy_stat*>** per_chr_per_bin_stats = new vector<t_per_bin_10x_discordancy_stat*>*[chr_ids->size() + 2];
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Generating bins on %s                   \r", chr_ids->at(i_chr));

		char cur_chr_multimapp_fp[1000];
		sprintf(cur_chr_multimapp_fp, "%s/%s.bin", multimapp_dir, chr_ids->at(i_chr));
		int l_mapp = 0;
		double* cur_chr_multi_mapp = load_normalized_multimappability_profile(cur_chr_multimapp_fp, l_mapp);

		per_chr_per_bin_stats[i_chr] = new vector<t_per_bin_10x_discordancy_stat*>();
		int cur_bin_start = 1;
		while(cur_bin_start + l_bin < chr_lengths->at(i_chr))
		{
			double total_multimapp = 0;
			for(int i = cur_bin_start; i <= cur_bin_start+l_bin-1; i++)
			{
				if(i < l_mapp)
				{
					total_multimapp += cur_chr_multi_mapp[i];
				}
			} // i loop.
			total_multimapp /= l_bin;

			t_per_bin_10x_discordancy_stat* cur_bin_stats = new t_per_bin_10x_discordancy_stat();
			cur_bin_stats->cis_discordant_linked_read_barcodes = new vector<double>();
			cur_bin_stats->trans_discordant_linked_read_barcodes = new vector<double>();
			cur_bin_stats->n_discordant_reads = 0;
			cur_bin_stats->cur_bin_start = cur_bin_start;
			cur_bin_stats->avg_multimapp = total_multimapp;

			per_chr_per_bin_stats[i_chr]->push_back(cur_bin_stats);
			cur_bin_start += l_bin;
		} // cur_bin_start

		delete [] cur_chr_multi_mapp;
	} // i_chr loop.
	fprintf(stderr, "\n");

	fprintf(stderr, "Processing %s\n", parsed_linked_reads_fp);
	FILE* f_parsed_linked_reads = open_f(parsed_linked_reads_fp, "r");
	
	// This is the set of reads stratified per chromosome for this linked read group.
	vector<t_10x_linked_read_info*>** cur_per_chr_group_reads = new vector<t_10x_linked_read_info*>*[chr_ids->size()];
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		cur_per_chr_group_reads[i_chr] = new vector<t_10x_linked_read_info*>();
	} // i_chr loop.

	int n_processed_groups = 0;
	while(1)
	{
		char* cur_line = getline(f_parsed_linked_reads);
		if(cur_line == NULL)
		{
			break;
		}

		n_processed_groups++;

		if(n_processed_groups % (100*1000) == 0)
		{
			fprintf(stderr, "Processing %d. GEM group.                   \r", n_processed_groups);
		}

		int i_cur_char = 0;
		char cur_tok[1000];

		// Read the barcode. Unused.
		t_string::get_next_token(cur_line, cur_tok, 1000, "\t", i_cur_char);
		double cur_barcode = atof(cur_tok);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
		fprintf(stderr, "New molecule: Barcode %.1f\n", cur_barcode);
}
		// Read the number of reads in this GEM group.
		t_string::get_next_token(cur_line, cur_tok, 1000, "\t", i_cur_char);
		int n_reads_per_cur_barcode = atoi(cur_tok);
		
		// Load the reads in this GEM group.
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
		fprintf(stderr, "Loading barcode: %.1f (%d)\n", cur_barcode, n_reads_per_cur_barcode);
}

		for(int read_i = 0; read_i < n_reads_per_cur_barcode; read_i++)
		{
			int cur_flag_self_chr_w_mate_chr = 0;
			int cur_posn = 0;
			int cur_mate_posn = 0;

			t_string::get_next_token(cur_line, cur_tok, 1000, "\t", i_cur_char);
			cur_flag_self_chr_w_mate_chr = atoi(cur_tok);

			t_string::get_next_token(cur_line, cur_tok, 1000, "\t", i_cur_char);
			cur_posn = atoi(cur_tok);

			t_string::get_next_token(cur_line, cur_tok, 1000, "\t", i_cur_char);
			cur_mate_posn = atoi(cur_tok);

			// Get this read's self chromosome index.
			int cur_read_chr_i = get_chr_i_per_group_entry(cur_flag_self_chr_w_mate_chr, true);
			int cur_read_flag = get_flag_per_group_entry(cur_flag_self_chr_w_mate_chr);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Read: CHR_i: %d, FLAG: %d (%d); POSNS: %d, %d\n",
					cur_read_chr_i,
					cur_read_flag, cur_flag_self_chr_w_mate_chr,
					cur_posn, cur_mate_posn);
			getc(stdin);
}

			t_10x_linked_read_info* cur_read_info = new t_10x_linked_read_info();
			cur_read_info->flag = cur_read_flag;
			cur_read_info->posn = cur_posn;
			cur_per_chr_group_reads[cur_read_chr_i]->push_back(cur_read_info);
		} // read_i loop.

		// Go over each chromosome and evaluate whether there is discordancy on any of them.		
		vector<int>* chr_i_per_chr_w_fragments = new vector<int>();
		int min_l_over_chromosomes = 1000*1000*100;
		int** per_chrom_fragment_start_ends = new int*[chr_ids->size() + 2];

		for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Processing %s: %d reads.\n", chr_ids->at(i_chr), (int)(cur_per_chr_group_reads[i_chr]->size()));
}

			per_chrom_fragment_start_ends[i_chr] = new int[2];
			per_chrom_fragment_start_ends[i_chr][0] = 0;
			per_chrom_fragment_start_ends[i_chr][1] = 0;

			// Make sure that we see good read support on this chromosome for this molecular fragment.
			if(cur_per_chr_group_reads[i_chr]->size() > min_n_supporting_reads_in_broken_fragment)
			{
				sort(cur_per_chr_group_reads[i_chr]->begin(), cur_per_chr_group_reads[i_chr]->end(), sort_10x_linked_read_info_per_position);

				// Set the start and end.
				per_chrom_fragment_start_ends[i_chr][0] = cur_per_chr_group_reads[i_chr]->at(0)->posn;
				per_chrom_fragment_start_ends[i_chr][1] = cur_per_chr_group_reads[i_chr]->back()->posn;

				// This is the length of the molecule on this current chromosome.
				int l_molecule_on_cur_chr = cur_per_chr_group_reads[i_chr]->back()->posn - cur_per_chr_group_reads[i_chr]->at(0)->posn + 1;

				// Update the minimum fragment length over the chromosomes.
				if(min_l_over_chromosomes > l_molecule_on_cur_chr)
				{
					min_l_over_chromosomes = l_molecule_on_cur_chr;
				}

				// If there is a good molecular length on this chromosome, add this.
				if(l_molecule_on_cur_chr > min_supporting_l_broken_fragment)
				{
					chr_i_per_chr_w_fragments->push_back(i_chr);
				}

				// Get the maximum interread distance on the molecule.
				vector<int>* interread_distances = new vector<int>();
				for(int i_read = 1;
					i_read < cur_per_chr_group_reads[i_chr]->size(); 
					i_read++)
				{
					int cur_dist = cur_per_chr_group_reads[i_chr]->at(i_read)->posn - cur_per_chr_group_reads[i_chr]->at(i_read - 1)->posn;
					if(cur_dist < 0)
					{
						fprintf(stderr, "Sanity check failed; inter-read distance negative for barcode %.1f: %d:%d, %d:%d\n", 
								cur_barcode, 
								i_chr, cur_per_chr_group_reads[i_chr]->at(i_read)->posn, 
								i_chr, cur_per_chr_group_reads[i_chr]->at(i_read-1)->posn);
						exit(0);
					}

					interread_distances->push_back(cur_dist);
				} // i_read loop.

				// Sort the distances.
				sort(interread_distances->begin(), interread_distances->end());

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				fprintf(stderr, "Barcode %.1f: max inter-read distance: %d (%d)\n", cur_barcode, 
					interread_distances->back(), 
					interread_distances->at(interread_distances->size()-2));
				getc(stdin);
}

				// We require at least 2 reads (we already have more), and significant discordancy between furthest reads.
				int n_distances = interread_distances->size();
				if(n_distances > 1 &&
					interread_distances->back() > min_l_discordant_interread_dist_on_molecule &&
					interread_distances->back() > min_mult_d12_per_distance * interread_distances->at(interread_distances->size()-2))
				{
					// Find the read that has the largest inter-read distance.
					for(int i_read = 1;
						i_read < cur_per_chr_group_reads[i_chr]->size(); 
						i_read++)
					{
						int cur_dist = cur_per_chr_group_reads[i_chr]->at(i_read)->posn - cur_per_chr_group_reads[i_chr]->at(i_read - 1)->posn;
						int n_read_2_left = i_read;
						int n_read_2_right = cur_per_chr_group_reads[i_chr]->size() - i_read + 1;

						// Make sure there is a good number of reads to left and to right.
						if(cur_dist == interread_distances->back() && 
							n_read_2_left > min_n_supporting_reads_in_broken_fragment &&	// Require at least 3 reads on each side.
							n_read_2_right > min_n_supporting_reads_in_broken_fragment)		// Require at least 3 reads on each side.
						{
							// This is the good breakpoint; add both sides to the block.
							int i_left_bin = (cur_per_chr_group_reads[i_chr]->at(i_read-1)->posn / l_bin);
							if(i_left_bin > per_chr_per_bin_stats[i_chr]->size()+10)
							{
								fprintf(stderr, "Sanity check failed: %d. bin requested out of %d bins.\n", 
									i_left_bin, per_chr_per_bin_stats[i_chr]->size());
								exit(0);
							}

							if(i_left_bin >= per_chr_per_bin_stats[i_chr]->size())
							{
								i_left_bin = per_chr_per_bin_stats[i_chr]->size() - 1;
							}

							// Update the molecular barcodes for this putative breakpoint.
							per_chr_per_bin_stats[i_chr]->at(i_left_bin)->cis_discordant_linked_read_barcodes->push_back(cur_barcode);

							// Add right junction information.
							int i_right_bin = (cur_per_chr_group_reads[i_chr]->at(i_read)->posn / l_bin);
							if(i_right_bin > per_chr_per_bin_stats[i_chr]->size()+10)
							{
								fprintf(stderr, "Sanity check failed: %d. bin requested out of %d bins.\n", 
									i_right_bin, per_chr_per_bin_stats[i_chr]->size());
								exit(0);
							}

							if(i_right_bin >= per_chr_per_bin_stats[i_chr]->size())
							{
								i_right_bin = per_chr_per_bin_stats[i_chr]->size() - 1;
							}

							// Update the molecular barcodes for this putative breakpoint.
							per_chr_per_bin_stats[i_chr]->at(i_right_bin)->cis_discordant_linked_read_barcodes->push_back(cur_barcode);
							break;
						} // i_read loop.
					} // i_read loop.
				} // discordancy check.

				delete interread_distances;
			} // multiple reads check.
		} // i_chr loop: cis BP check.

		// Check fragments on multiple chromosomes.
		if(chr_i_per_chr_w_fragments->size() > 1 &&						// Require multiple chromosomes.
			min_l_over_chromosomes > min_supporting_l_broken_fragment)	// Require at least 300 bp span on each chromosome.
		{
			// Discordancy by fragments on multiple chromosomes.
			for(int i_chr_i = 0; i_chr_i < chr_i_per_chr_w_fragments->size(); i_chr_i++)
			{
				int i_chr = chr_i_per_chr_w_fragments->at(i_chr_i);

				sort(cur_per_chr_group_reads[i_chr]->begin(), cur_per_chr_group_reads[i_chr]->end(), sort_10x_linked_read_info_per_position);

				// Update the counts at the beginning and end of this fragment.
				int i_left_bin = (cur_per_chr_group_reads[i_chr]->at(0)->posn / l_bin);
				if(i_left_bin >= per_chr_per_bin_stats[i_chr]->size())
				{
					i_left_bin = per_chr_per_bin_stats[i_chr]->size() - 1;
				}

				// Update the molecular barcodes for this putative breakpoint.
				per_chr_per_bin_stats[i_chr]->at(i_left_bin)->trans_discordant_linked_read_barcodes->push_back(cur_barcode);

				// Add right junction information.
				int i_right_bin = (cur_per_chr_group_reads[i_chr]->back()->posn / l_bin);
				if(i_right_bin >= per_chr_per_bin_stats[i_chr]->size())
				{
					i_right_bin = per_chr_per_bin_stats[i_chr]->size() - 1;
				}

				// Update the molecular barcodes for this putative breakpoint.
				per_chr_per_bin_stats[i_chr]->at(i_right_bin)->trans_discordant_linked_read_barcodes->push_back(cur_barcode);
				break;
			} // i_chr_i loop.
		} // trans BP check.

		// Free the read memory.
		delete chr_i_per_chr_w_fragments;
		for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			// Free read memory.
			for(int i_read = 0; i_read < cur_per_chr_group_reads[i_chr]->size(); i_read++)
			{
				delete cur_per_chr_group_reads[i_chr]->at(i_read);
			} // i_read loop.

			// Clear the list of reads.
			cur_per_chr_group_reads[i_chr]->clear();

			delete [] per_chrom_fragment_start_ends[i_chr];
		} // i_chr loop.
		delete [] per_chrom_fragment_start_ends;

		// Free this line.
		delete [] cur_line;
	} // file reading loop.
	fclose(f_parsed_linked_reads);

	fprintf(stderr, "Saving the binned linked read stats to %s\n", op_fp);
	FILE* f_op = open_f(op_fp, "w");
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Saving information for %s\n", chr_ids->at(i_chr));
		for(int i_bin = 0; i_bin < per_chr_per_bin_stats[i_chr]->size(); i_bin++)
		{
			if(per_chr_per_bin_stats[i_chr]->at(i_bin)->cis_discordant_linked_read_barcodes->size() > 0 ||
				per_chr_per_bin_stats[i_chr]->at(i_bin)->trans_discordant_linked_read_barcodes->size() > 0)
			{
				fprintf(f_op, "%s\t%d\t%.4f\t%d\t%d", 
						chr_ids->at(i_chr), 
						per_chr_per_bin_stats[i_chr]->at(i_bin)->cur_bin_start, 
						per_chr_per_bin_stats[i_chr]->at(i_bin)->avg_multimapp, 
						(int)(per_chr_per_bin_stats[i_chr]->at(i_bin)->cis_discordant_linked_read_barcodes->size()),
						(int)(per_chr_per_bin_stats[i_chr]->at(i_bin)->trans_discordant_linked_read_barcodes->size()));

				for(int i_cis = 0; i_cis < per_chr_per_bin_stats[i_chr]->at(i_bin)->cis_discordant_linked_read_barcodes->size(); i_cis++)
				{
					fprintf(f_op, "\t%.1f", per_chr_per_bin_stats[i_chr]->at(i_bin)->cis_discordant_linked_read_barcodes->at(i_cis));
				} // i_cis loop.

				for(int i_trans = 0; i_trans < per_chr_per_bin_stats[i_chr]->at(i_bin)->trans_discordant_linked_read_barcodes->size(); i_trans++)
				{
					fprintf(f_op, "\t%.1f", per_chr_per_bin_stats[i_chr]->at(i_bin)->trans_discordant_linked_read_barcodes->at(i_trans));
				} // i_cis loop.

				fprintf(f_op, "\n");
			} // check if there are cis or trans discordant reads.
		} // i_bin loop.
	} // i_chr loop.
	fclose(f_op);
} // compute_discordancy_stats_per_10x_linked_reads

void sort_PE_reads_file_per_id_in_memory(char* mapped_reads_fp,
	char* sorted_op_fp)
{
	vector<char*>* read_lines = buffer_file(mapped_reads_fp);
	if(read_lines == NULL)
	{
		fprintf(stderr, "Could not open %s\n", mapped_reads_fp);
		exit(0);
	}

	sort_read_lines_per_id_in_place(read_lines);
		
	FILE* f_cur_chr_op = open_f(sorted_op_fp, "w");
	for(int i_read = 0; i_read < (int)read_lines->size(); i_read++)
	{
		fprintf(f_cur_chr_op, "%s\n", read_lines->at(i_read));
	}  // i_read loop.
	fclose(f_cur_chr_op);
} // sort_PE_reads_file_per_id_in_memory

/*
Sort the reads per id per chromosome first, then do external sorting on all the files.
*/
void label_merge_external_sort_PE_reads_per_id(char* first_mapped_reads_dir,
	char* last_mapped_reads_dir,
	char* temp_sorted_op_dir,
	char* sorted_op_fp)
{
	// Load all the first pair id's.
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", first_mapped_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	if(chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_fp);
		exit(0);
	}
	
	vector<FILE*>* per_chromosome_fps = new vector<FILE*>();
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Load and process the first pairs.
		char cur_chr_reads_fp[1000];
		//sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", first_mapped_reads_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(first_mapped_reads_dir, chr_ids->at(i_chr), cur_chr_reads_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", first_mapped_reads_dir, chr_ids->at(i_chr));
			exit(0);
		}

		vector<char*>* cur_chr_first_read_lines = buffer_file(cur_chr_reads_fp);
		if(cur_chr_first_read_lines == NULL)
		{
			fprintf(stderr, "Could not open %s\n", cur_chr_reads_fp);
			exit(0);
		}

		// Add the end id to the line.
		fprintf(stderr, "Copying the mate id to first pair reads.\n");
		for(int i_l = 0; i_l < (int)cur_chr_first_read_lines->size(); i_l++)
		{
			char* cur_line_w_pair_id = new char[t_string::string_length(cur_chr_first_read_lines->at(i_l)) + 4];
			sprintf(cur_line_w_pair_id, "%s 0", cur_chr_first_read_lines->at(i_l));
			delete [] cur_chr_first_read_lines->at(i_l);

			// Replace the last read line with the line with id.
			cur_chr_first_read_lines->at(i_l) = cur_line_w_pair_id;
		} // i_l loop.

		// Load and process the last pairs.
		//sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", last_mapped_reads_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(last_mapped_reads_dir, chr_ids->at(i_chr), cur_chr_reads_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", last_mapped_reads_dir, chr_ids->at(i_chr));
			exit(0);
		}

		vector<char*>* cur_chr_last_read_lines = buffer_file(cur_chr_reads_fp);

		// Add the end id to the line.
		fprintf(stderr, "Copying the mate id to last pair reads.\n");
		for(int i_l = 0; i_l < (int)cur_chr_last_read_lines->size(); i_l++)
		{
			char* cur_line_w_pair_id = new char[t_string::string_length(cur_chr_last_read_lines->at(i_l)) + 4];
			sprintf(cur_line_w_pair_id, "%s 1", cur_chr_last_read_lines->at(i_l));
			delete [] cur_chr_last_read_lines->at(i_l);

			// Replace the last read line with the line with id.
			cur_chr_last_read_lines->at(i_l) = cur_line_w_pair_id;
		} // i_l loop.

		// Add the last reads to the first reads.
		cur_chr_first_read_lines->insert(cur_chr_first_read_lines->end(), 
			cur_chr_last_read_lines->begin(), cur_chr_last_read_lines->end());

		// We do not need this vector any more.
		delete cur_chr_last_read_lines;

		fprintf(stderr, "Sorting %ld reads per id\n", cur_chr_first_read_lines->size());
		sort_read_lines_per_id_in_place(cur_chr_first_read_lines);
		
		char cur_chr_op_fp[1000];
		//sprintf(cur_chr_op_fp, "%s/%s_mapped_reads.txt", temp_sorted_op_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(temp_sorted_op_dir, chr_ids->at(i_chr), cur_chr_op_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", temp_sorted_op_dir, chr_ids->at(i_chr));
			exit(0);
		}

		fprintf(stderr, "Sorted %s: %ld reads. Dumping to %s.\n", chr_ids->at(i_chr), cur_chr_first_read_lines->size(), cur_chr_op_fp);		
		FILE* f_cur_chr_op = open_f(cur_chr_op_fp, "w");
		for(int i_read = 0; i_read < (int)cur_chr_first_read_lines->size(); i_read++)
		{
			fprintf(f_cur_chr_op, "%s\t%s\n", cur_chr_first_read_lines->at(i_read), chr_ids->at(i_chr));
		}  // i_read loop.
		fclose(f_cur_chr_op);

		//delete cur_chr_first_read_lines_pair_id;
		t_string::clean_string_list(cur_chr_first_read_lines);

		// Add the current chromosome.
		per_chromosome_fps->push_back(open_f(cur_chr_op_fp, "r"));
	} // i_chr loop.

	fprintf(stderr, "Doing external sorting on %ld files.\n", per_chromosome_fps->size());

	// Now open the per chromosome sorted reads and start doing external sort.
	vector<char*>* cur_ids_per_chrom = new vector<char*>();
	vector<char*>* cur_lines_per_chrom = new vector<char*>();
	int i_cur_first_id = (int)per_chromosome_fps->size();
	char cur_read_id[1000];
	for(int i_f = 0; i_f < (int)per_chromosome_fps->size(); i_f++)
	{	
		char* cur_line = getline(per_chromosome_fps->at(i_f));
		sscanf(cur_line, "%s", cur_read_id);
		cur_lines_per_chrom->push_back(cur_line);
		cur_ids_per_chrom->push_back(t_string::copy_me_str(cur_read_id));
		
		if(i_cur_first_id == (int)per_chromosome_fps->size())
		{
			i_cur_first_id = i_f;
		}
		else
		{
			if(!t_string::sort_strings(cur_ids_per_chrom->at(i_cur_first_id), cur_ids_per_chrom->at(i_f)))
			{
				i_cur_first_id = i_f;
			}
		}
	} // i_f loop.

	// Dump the sorted reads.
	FILE* f_sorted = open_f(sorted_op_fp, "w");
	while(per_chromosome_fps->size() > 0)
	{
		// Dump the first, replace it find the new first.
		fprintf(f_sorted, "%s\n", cur_lines_per_chrom->at(i_cur_first_id));
		delete [] cur_lines_per_chrom->at(i_cur_first_id);
		delete [] cur_ids_per_chrom->at(i_cur_first_id);

		// Get the next line and copy the current id and the current line.
		cur_lines_per_chrom->at(i_cur_first_id) = getline(per_chromosome_fps->at(i_cur_first_id));

		if(cur_lines_per_chrom->at(i_cur_first_id) == NULL)
		{
			fprintf(stderr, "Finished %s\n", chr_ids->at(i_cur_first_id));

			// Close the file.
			fclose(per_chromosome_fps->at(i_cur_first_id));

			per_chromosome_fps->erase(per_chromosome_fps->begin() + i_cur_first_id);
			cur_lines_per_chrom->erase(cur_lines_per_chrom->begin() + i_cur_first_id);
			cur_ids_per_chrom->erase(cur_ids_per_chrom->begin() + i_cur_first_id);
			chr_ids->erase(chr_ids->begin() + i_cur_first_id);

			i_cur_first_id = (int)per_chromosome_fps->size();
		}
		else
		{
			sscanf(cur_lines_per_chrom->at(i_cur_first_id), "%s", cur_read_id);
			cur_ids_per_chrom->at(i_cur_first_id) = t_string::copy_me_str(cur_read_id);
		}

		// Update the current first.
		for(int i_f = 0; i_f < (int)per_chromosome_fps->size(); i_f++)
		{
			if(i_cur_first_id == (int)per_chromosome_fps->size())
			{
				i_cur_first_id = i_f;
			}
			else if(!t_string::sort_strings(cur_ids_per_chrom->at(i_cur_first_id), cur_ids_per_chrom->at(i_f)))
			{
				i_cur_first_id = i_f;
			}
		} // i_f loop.
	} // file reading loop.
	fclose(f_sorted);
} // label_merge_external_sort_PE_reads_per_id

//void label_PE_reads_file(char* mapped_reads_dir, char side_char, char* op_dir)
void label_PE_reads_file(char* reads_fp, char side_char, char* chr_id, char* op_fp)
{
	//// Load all the first pair id's.
	//char chr_ids_fp[1000];
	//sprintf(chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	//vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	//if(chr_ids == NULL)
	//{
	//	fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_fp);
	//	exit(0);
	//}
	
	//for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	//{
		fprintf(stderr, "Labeling reads on %s\n", reads_fp);
		//char cur_chr_reads_fp[1000];
		//sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, chr_ids->at(i_chr));
		FILE* f_cur_reads = open_f(reads_fp, "r");

		//char cur_chr_labeled_reads_fp[1000];
		//sprintf(cur_chr_labeled_reads_fp, "%s/%s_mapped_reads.txt", op_dir, chr_ids->at(i_chr));
		FILE* f_cur_labeled_reads = open_f(op_fp, "w");
		while(1)
		{
			char* cur_line = getline(f_cur_reads);
			if(cur_line == NULL)
			{
				break;
			}

			fprintf(f_cur_labeled_reads, "%s %c %s\n", cur_line, side_char, chr_id);

			delete [] cur_line;
		} // file readingl loop.
		fclose(f_cur_reads);
		fclose(f_cur_labeled_reads);
	//} // i_chr loop.
}

void external_sort_PE_reads_per_file_list(char* read_list_fp, char* sorted_op_fp)
{
	vector<char*>* read_fps = buffer_file(read_list_fp);

	vector<FILE*>* read_files = new vector<FILE*>();
	for(int i_f = 0; i_f < (int)read_fps->size(); i_f++)
	{
		// Add the current chromosome.
		read_files->push_back(open_f(read_fps->at(i_f), "r"));
	} // i_chr loop.

	fprintf(stderr, "Doing external sorting on %ld files.\n", read_files->size());

	// Now open the per chromosome sorted reads and start doing external sort.
	vector<char*>* cur_ids_per_chrom = new vector<char*>();
	vector<char*>* cur_lines_per_chrom = new vector<char*>();
	int i_cur_first_id = (int)read_files->size();
	char cur_read_id[1000];
	for(int i_f = 0; i_f < (int)read_files->size(); i_f++)
	{	
		char* cur_line = getline(read_files->at(i_f));
		sscanf(cur_line, "%s", cur_read_id);
		cur_lines_per_chrom->push_back(cur_line);
		cur_ids_per_chrom->push_back(t_string::copy_me_str(cur_read_id));
		
		if(i_cur_first_id == (int)read_files->size())
		{
			i_cur_first_id = i_f;
		}
		else
		{
			if(!t_string::sort_strings(cur_ids_per_chrom->at(i_cur_first_id), cur_ids_per_chrom->at(i_f)))
			{
				i_cur_first_id = i_f;
			}
		}
	} // i_f loop.

	// Dump the sorted reads.
	FILE* f_sorted = open_f(sorted_op_fp, "w");
	while(read_files->size() > 0)
	{
		// Dump the first, replace it find the new first.
		fprintf(f_sorted, "%s\n", cur_lines_per_chrom->at(i_cur_first_id));
		delete [] cur_lines_per_chrom->at(i_cur_first_id);
		delete [] cur_ids_per_chrom->at(i_cur_first_id);

		// Get the next line and copy the current id and the current line.
		cur_lines_per_chrom->at(i_cur_first_id) = getline(read_files->at(i_cur_first_id));

		if(cur_lines_per_chrom->at(i_cur_first_id) == NULL)
		{
			fprintf(stderr, "Finished %s\n", read_fps->at(i_cur_first_id));

			// Close the file.
			fclose(read_files->at(i_cur_first_id));

			read_files->erase(read_files->begin() + i_cur_first_id);
			cur_lines_per_chrom->erase(cur_lines_per_chrom->begin() + i_cur_first_id);
			cur_ids_per_chrom->erase(cur_ids_per_chrom->begin() + i_cur_first_id);
			read_fps->erase(read_fps->begin() + i_cur_first_id);

			i_cur_first_id = (int)read_files->size();
		}
		else
		{
			sscanf(cur_lines_per_chrom->at(i_cur_first_id), "%s", cur_read_id);
			cur_ids_per_chrom->at(i_cur_first_id) = t_string::copy_me_str(cur_read_id);
		}

		// Update the current first.
		for(int i_f = 0; i_f < (int)read_files->size(); i_f++)
		{
			if(i_cur_first_id == (int)read_files->size())
			{
				i_cur_first_id = i_f;
			}
			else if(!t_string::sort_strings(cur_ids_per_chrom->at(i_cur_first_id), cur_ids_per_chrom->at(i_f)))
			{
				i_cur_first_id = i_f;
			}
		} // i_f loop.
	} // file reading loop.
	fclose(f_sorted);

	t_string::clean_string_list(read_fps);
}

void label_merge_external_sort_SE_reads_per_id(char* mapped_reads_dir,
	char* temp_sorted_op_dir,
	char* sorted_op_fp)
{
	// Load all the first pair id's.
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	if(chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_fp);
		exit(0);
	}
	
	vector<FILE*>* per_chromosome_fps = new vector<FILE*>();
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Load and process the first pairs.
		char cur_chr_reads_fp[1000];
		//sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(mapped_reads_dir, chr_ids->at(i_chr), cur_chr_reads_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", temp_sorted_op_dir, chr_ids->at(i_chr));
			exit(0);
		}

		vector<char*>* cur_chr_first_read_lines = buffer_file(cur_chr_reads_fp);
		if(cur_chr_first_read_lines == NULL)
		{
			fprintf(stderr, "Could not open %s\n", cur_chr_reads_fp);
			exit(0);
		}

		// Add the end id to the line.
		fprintf(stderr, "Copying the mate id to first pair reads.\n");
		vector<char*>* cur_chr_read_lines_pair_id = new vector<char*>();
		for(int i_l = 0; i_l < (int)cur_chr_first_read_lines->size(); i_l++)
		{
			char* cur_line_w_pair_id = new char[t_string::string_length(cur_chr_first_read_lines->at(i_l)) + 4];
			sprintf(cur_line_w_pair_id, "%s 0", cur_chr_first_read_lines->at(i_l));
			cur_chr_read_lines_pair_id->push_back(cur_line_w_pair_id);
			delete [] cur_chr_first_read_lines->at(i_l);
		} // i_l loop.
		delete cur_chr_first_read_lines;

		fprintf(stderr, "Sorting %ld reads per id\n", cur_chr_read_lines_pair_id->size());
		vector<char*>* sorted_read_lines = sort_read_lines_per_id(cur_chr_read_lines_pair_id);
		
		char cur_chr_op_fp[1000];
		//sprintf(cur_chr_op_fp, "%s/%s_mapped_reads.txt", temp_sorted_op_dir, chr_ids->at(i_chr));
		if (!set_preprocessed_read_file_path_per_dir_chr(temp_sorted_op_dir, chr_ids->at(i_chr), cur_chr_op_fp))
		{
			fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", temp_sorted_op_dir, chr_ids->at(i_chr));
			exit(0);
		}

		fprintf(stderr, "Sorted %s: %ld reads. Dumping to %s.\n", chr_ids->at(i_chr), sorted_read_lines->size(), cur_chr_op_fp);
		FILE* f_cur_chr_op = open_f(cur_chr_op_fp, "w");
		for(int i_read = 0; i_read < (int)sorted_read_lines->size(); i_read++)
		{
			fprintf(f_cur_chr_op, "%s\t%s\n", sorted_read_lines->at(i_read), chr_ids->at(i_chr));
		}  // i_read loop.
		fclose(f_cur_chr_op);

		delete cur_chr_read_lines_pair_id;

		// Clean the read lines.
		t_string::clean_string_list(sorted_read_lines);

		// Add the current chromosome.
		per_chromosome_fps->push_back(open_f(cur_chr_op_fp, "r"));
	} // i_chr loop.

	fprintf(stderr, "Doing external sorting on %ld files.\n", per_chromosome_fps->size());

	// Now open the per chromosome sorted reads and start doing external sort.
	vector<char*>* cur_ids_per_chrom = new vector<char*>();
	vector<char*>* cur_lines_per_chrom = new vector<char*>();
	int i_cur_first_id = (int)per_chromosome_fps->size();
	char cur_read_id[1000];
	for(int i_f = 0; i_f < (int)per_chromosome_fps->size(); i_f++)
	{	
		char* cur_line = getline(per_chromosome_fps->at(i_f));
		sscanf(cur_line, "%s", cur_read_id);
		cur_lines_per_chrom->push_back(cur_line);
		cur_ids_per_chrom->push_back(t_string::copy_me_str(cur_read_id));
		
		if(i_cur_first_id == (int)per_chromosome_fps->size())
		{
			i_cur_first_id = i_f;
		}
		else
		{
			if(!t_string::sort_strings(cur_ids_per_chrom->at(i_cur_first_id), cur_ids_per_chrom->at(i_f)))
			{
				i_cur_first_id = i_f;
			}
		}
	} // i_f loop.

	// Dump the sorted reads.
	FILE* f_sorted = open_f(sorted_op_fp, "w");
	while(per_chromosome_fps->size() > 0)
	{
		// Dump the first, replace it find the new first.
		fprintf(f_sorted, "%s\n", cur_lines_per_chrom->at(i_cur_first_id));
		delete [] cur_lines_per_chrom->at(i_cur_first_id);
		delete [] cur_ids_per_chrom->at(i_cur_first_id);

		// Get the next line and copy the current id and the current line.
		cur_lines_per_chrom->at(i_cur_first_id) = getline(per_chromosome_fps->at(i_cur_first_id));

		if(cur_lines_per_chrom->at(i_cur_first_id) == NULL)
		{
			fprintf(stderr, "Finished %s\n", chr_ids->at(i_cur_first_id));

			// Close the file.
			fclose(per_chromosome_fps->at(i_cur_first_id));

			per_chromosome_fps->erase(per_chromosome_fps->begin() + i_cur_first_id);
			cur_lines_per_chrom->erase(cur_lines_per_chrom->begin() + i_cur_first_id);
			cur_ids_per_chrom->erase(cur_ids_per_chrom->begin() + i_cur_first_id);
			chr_ids->erase(chr_ids->begin() + i_cur_first_id);

			i_cur_first_id = (int)per_chromosome_fps->size();
		}
		else
		{
			sscanf(cur_lines_per_chrom->at(i_cur_first_id), "%s", cur_read_id);
			cur_ids_per_chrom->at(i_cur_first_id) = t_string::copy_me_str(cur_read_id);
		}

		// Update the current first.
		for(int i_f = 0; i_f < (int)per_chromosome_fps->size(); i_f++)
		{
			if(i_cur_first_id == (int)per_chromosome_fps->size())
			{
				i_cur_first_id = i_f;
			}
			else if(!t_string::sort_strings(cur_ids_per_chrom->at(i_cur_first_id), cur_ids_per_chrom->at(i_f)))
			{
				i_cur_first_id = i_f;
			}
		} // i_f loop.
	} // file reading loop.
	fclose(f_sorted);
}

#define L_CHROM (250*1000*1000)
//void generate_discordant_PE_read_connections_per_sorted_PE_reads_file(char* chr_ids_fp, char* sorted_pe_reads_fp, 
//	int l_bin, 
//	int min_l_same_chr_separation, int max_concordant_mapping_separation,
//	char valid_first_mapper_strand,
//	char valid_last_mapper_strand)
//{
//	// Allocate the bins for counting the connections.
//	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
//	int n_bins_per_chr = (L_CHROM / l_bin) + 10;
//
//	// Allocate the bins in each chromosome.
//	vector<vector<t_annot_region*>*>* bin_regions_per_chr = new vector<vector<t_annot_region*>*>();
//	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
//	{		
//		fprintf(stderr, "Allocating bins in %s\n", chr_ids->at(i_chr));
//		vector<t_annot_region*>* bin_regions_per_cur_chr = new vector<t_annot_region*>();
//
//		int cur_bin_start = 1;
//		for(int i_bin = 0; i_bin < n_bins_per_chr; i_bin++)
//		{
//			t_annot_region* cur_bin_region = get_empty_region();
//			cur_bin_region->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
//			cur_bin_region->start = cur_bin_start;
//			cur_bin_region->end = cur_bin_start + l_bin;
//			vector<t_connection_info_entry*>* cur_bin_connecting_bins = new vector<t_connection_info_entry*>();
//			cur_bin_region->data = cur_bin_connecting_bins;
//
//			bin_regions_per_cur_chr->push_back(cur_bin_region);
//
//			cur_bin_start += l_bin;
//		} // i_bin loop.
//
//		fprintf(stderr, "%ld regions (%d)\n", bin_regions_per_cur_chr->size(), cur_bin_start);
//		bin_regions_per_chr->push_back(bin_regions_per_cur_chr);
//	} // i_chr loop.
//
//	FILE* f_sorted_pe_reads = open_f(sorted_pe_reads_fp, "r");
//	char prev_read_id[1000];
//	prev_read_id[0] = 0;
//
//	vector<char*>* cur_first_mapper_read_ids = new vector<char*>();
//	vector<int>* cur_first_mapper_indices = new vector<int>();
//	vector<char*>* cur_first_mapper_chrs = new vector<char*>();
//	vector<char>* cur_first_mapper_strands = new vector<char>();
//	vector<char*>* cur_last_mapper_read_ids = new vector<char*>();
//	vector<int>* cur_last_mapper_indices = new vector<int>();
//	vector<char*>* cur_last_mapper_chrs = new vector<char*>();
//	vector<char>* cur_last_mapper_strands = new vector<char>();
//
//	// Get the current list of connections: Get the consecutive lines with same read id.
//	// MARILYN_0005:1:2:1272:993#0 76M R 102495555 0   5
//	int n_lines_processed = 0;
//	int n_concordant = 0;
//	int n_discordant = 0;
//	while(1)
//	{
//		if(n_lines_processed % 1000000 == 0)
//		{
//			fprintf(stderr, "Processed %d (%d, %d) read.           \r", n_lines_processed, n_concordant, n_discordant);
//		}
//
//		// Read the current read line.
//		char* cur_read_line = getline(f_sorted_pe_reads);
//		if(cur_read_line == NULL)
//		{
//			break;
//		}
//
//		n_lines_processed++;
//
//		char cur_read_id[1000];
//		char cur_mapping_str[1000];
//		char cur_strand;
//		int cur_start;
//		int cur_pair_index;
//		char cur_chr_id[100];
//
//		if(sscanf(cur_read_line, "%s %s %c %d %d %s", 
//			cur_read_id, cur_mapping_str, &cur_strand, &cur_start, &cur_pair_index, cur_chr_id) != 6)
//		{
//			fprintf(stderr, "Could not parse: %s\n", cur_read_line);
//			exit(0);
//		}
//
//		// Does this read have the same id as previous ones?
//		if(t_string::compare_strings(prev_read_id, cur_read_id))
//		{
//			if(cur_pair_index == 0)
//			{
//				cur_first_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
//				cur_first_mapper_indices->push_back(cur_start);
//				cur_first_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
//				cur_first_mapper_strands->push_back(cur_strand);
//			}
//			else
//			{
//				cur_last_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
//				cur_last_mapper_indices->push_back(cur_start);
//				cur_last_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
//				cur_last_mapper_strands->push_back(cur_strand);
//			}
//		}
//		else
//		{
//			bool found_concordant_mapping = false;
//			for(int i_first = 0; 
//				!found_concordant_mapping && i_first < (int)cur_first_mapper_indices->size(); 
//				i_first++)
//			{
//				for(int i_last = 0; 
//					!found_concordant_mapping && i_last < (int)cur_last_mapper_indices->size(); 
//					i_last++)
//				{
//					// Update the current pair.
//					char first_mapper_strand = cur_first_mapper_strands->at(i_first);
//					char last_mapper_strand = cur_last_mapper_strands->at(i_last);
//					char* first_mapper_chr = cur_first_mapper_chrs->at(i_first);
//					char* last_mapper_chr = cur_last_mapper_chrs->at(i_last);
//
//					// Is this the concordant mapping for this current read?
//					if(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
//						abs(cur_first_mapper_indices->at(i_first) - cur_last_mapper_indices->at(i_last)) < max_concordant_mapping_separation &&
//						((first_mapper_strand == valid_first_mapper_strand && last_mapper_strand == valid_last_mapper_strand) ||
//						(first_mapper_strand == valid_last_mapper_strand && last_mapper_strand == valid_first_mapper_strand)))
//					{
//						found_concordant_mapping = true;
//					}
//				} // i_last loop.
//			} // i_first loop.
//
//			// Process the curent list of entries if there was no concordant mapping.
//			if(!found_concordant_mapping)
//			{
//				n_discordant++;
//
//if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
//{
//                if(cur_first_mapper_indices->size() > 0 && cur_last_mapper_indices->size() > 0)
//                {
//					fprintf(stderr, "First reads:\n");
//					for(int i_first = 0; 
//						i_first < (int)cur_first_mapper_indices->size(); 
//						i_first++)
//					{
//						fprintf(stderr, "%s: %s: %d (%c)\n", cur_first_mapper_read_ids->at(i_first), 
//							cur_first_mapper_chrs->at(i_first), 
//							cur_first_mapper_indices->at(i_first), 
//							cur_first_mapper_strands->at(i_first));
//					} // i_first loop.
//
//					fprintf(stderr, "Last reads:\n");
//					for(int i_last = 0; 
//						i_last < (int)cur_last_mapper_indices->size(); 
//						i_last++)
//					{
//						fprintf(stderr, "%s: %s: %d (%c)\n", cur_last_mapper_read_ids->at(i_last), 
//							cur_last_mapper_chrs->at(i_last),
//							cur_last_mapper_indices->at(i_last), 
//							cur_last_mapper_strands->at(i_last));
//					} // i_last loop.
//					getc(stdin);
//				}
//} // check for dumping reads.
//
//				for(int i_first = 0; 
//					i_first < (int)cur_first_mapper_indices->size(); 
//					i_first++)
//				{
//					for(int i_last = 0; 
//						i_last < (int)cur_last_mapper_indices->size(); 
//						i_last++)
//					{
//						// Update the current pair.
//						char first_mapper_strand = cur_first_mapper_strands->at(i_first);
//						char last_mapper_strand = cur_last_mapper_strands->at(i_last);
//
//						char* first_mapper_chr = cur_first_mapper_chrs->at(i_first);
//						char* last_mapper_chr = cur_last_mapper_chrs->at(i_last);
//
//						// Since we know that the pairs are discordant, we do not need any more filters.
//						//if(first_mapper_strand != last_mapper_strand &&			
//						//	(!t_string::compare_strings(first_mapper_chr, last_mapper_chr) || 
//						//	(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
//						//	abs(cur_first_mapper_indices->at(i_first) - cur_last_mapper_indices->at(i_last)) > min_l_same_chr_separation)))
//						//if(!t_string::compare_strings(first_mapper_chr, last_mapper_chr))
//						{
//							// Get the chr indices for the mapping target bins.
//							int i_first_chr_id = t_string::get_i_str(chr_ids, cur_first_mapper_chrs->at(i_first));
//							int i_last_chr_id = t_string::get_i_str(chr_ids, cur_last_mapper_chrs->at(i_last));
//
//							t_annot_region* cur_first_bin = bin_regions_per_chr->at(i_first_chr_id)->at(cur_first_mapper_indices->at(i_first) / l_bin);
//							t_annot_region* cur_last_bin = bin_regions_per_chr->at(i_last_chr_id)->at(cur_last_mapper_indices->at(i_last) / l_bin);
//
//							vector<t_connection_info_entry*>* cur_first_bin_bins_list = (vector<t_connection_info_entry*>*)(cur_first_bin->data);
//							vector<t_connection_info_entry*>* cur_last_bin_bins_list = (vector<t_connection_info_entry*>*)(cur_last_bin->data);
//
//							// Add the regions to the connection lists, if they do not exist, add them, otherwise include them.
//							bool found_connection1 = false;
//							for(int i_conn1 = 0; 
//								i_conn1 < cur_first_bin_bins_list->size() &&
//								!found_connection1; 
//								i_conn1++)
//							{
//								if(t_string::compare_strings(cur_last_bin->chrom, cur_first_bin_bins_list->at(i_conn1)->connecting_region->chrom) &&
//									cur_last_bin->start == cur_first_bin_bins_list->at(i_conn1)->connecting_region->start &&
//									cur_last_bin->end == cur_first_bin_bins_list->at(i_conn1)->connecting_region->end)
//								{
//									cur_first_bin_bins_list->at(i_conn1)->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//									cur_first_bin_bins_list->at(i_conn1)->first_read_strands->push_back(first_mapper_strand);
//									cur_first_bin_bins_list->at(i_conn1)->last_read_strands->push_back(last_mapper_strand);
//									found_connection1 = true;
//								}
//							} // i_conn1 loop.
//
//							if(!found_connection1)
//							{
//								t_connection_info_entry* new_first_connection_entry = new t_connection_info_entry();
//								new_first_connection_entry->connecting_region = cur_last_bin;
//								new_first_connection_entry->read_ids = new vector<char*>();
//								new_first_connection_entry->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//								new_first_connection_entry->first_read_strands = new vector<char>();
//								new_first_connection_entry->first_read_strands->push_back(first_mapper_strand);
//								new_first_connection_entry->last_read_strands = new vector<char>();
//								new_first_connection_entry->last_read_strands->push_back(last_mapper_strand);
//								cur_first_bin_bins_list->push_back(new_first_connection_entry);
//							}
//
//							// Look for the first bin in the connection list of last bin.
//							// If the bins are the same bin, do not run following as it is redundant.
//							bool found_connection2 = false;
//							if(cur_first_bin != cur_last_bin)
//							{
//								for(int i_conn1 = 0; 
//									i_conn1 < cur_last_bin_bins_list->size() &&
//									!found_connection2; 
//									i_conn1++)
//								{
//									if(t_string::compare_strings(cur_first_bin->chrom, cur_last_bin_bins_list->at(i_conn1)->connecting_region->chrom) &&
//										cur_first_bin->start ==  cur_last_bin_bins_list->at(i_conn1)->connecting_region->start &&
//										cur_first_bin->end == cur_last_bin_bins_list->at(i_conn1)->connecting_region->end)
//									{
//										cur_last_bin_bins_list->at(i_conn1)->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//										cur_last_bin_bins_list->at(i_conn1)->first_read_strands->push_back(first_mapper_strand);
//										cur_last_bin_bins_list->at(i_conn1)->last_read_strands->push_back(last_mapper_strand);
//										found_connection2 = true;
//									}
//								} // i_conn1 loop.
//
//								if(!found_connection2)
//								{
//									t_connection_info_entry* new_last_connection_entry = new t_connection_info_entry();							
//									new_last_connection_entry->connecting_region = cur_first_bin;
//									new_last_connection_entry->read_ids = new vector<char*>();
//									new_last_connection_entry->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//									new_last_connection_entry->first_read_strands = new vector<char>();
//									new_last_connection_entry->first_read_strands->push_back(first_mapper_strand);
//									new_last_connection_entry->last_read_strands = new vector<char>();
//									new_last_connection_entry->last_read_strands->push_back(last_mapper_strand);
//									cur_last_bin_bins_list->push_back(new_last_connection_entry);
//								}
//							} // self check for cutting this loop.
//
//							if(cur_first_bin != cur_last_bin &&
//								found_connection2 != found_connection1)
//							{
//								fprintf(stderr, "The connections are not consistent: %s\t%d\t%d, %s\t%d\t%d\nLeft mapping: %d(%c), Right mapping: %d(%c)", 
//									cur_first_bin->chrom, cur_first_bin->start, cur_first_bin->end, 
//									cur_last_bin->chrom, cur_last_bin->start, cur_last_bin->end,
//									cur_first_mapper_indices->at(i_first), first_mapper_strand, 
//									cur_last_mapper_indices->at(i_last), last_mapper_strand);
//								exit(0);
//							}
//						} // extra discordance checks.
//					} // i_last loop.
//				} // i_first loop.
//			} // concordant read check.
//			else
//			{
//				n_concordant++;
//			}
//
//			// Process the current entry, clear the list, then add the entry for the current read.
//			cur_first_mapper_indices->clear();
//			cur_first_mapper_strands->clear();
//			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_first_mapper_chrs->size(); i_mapper_chr++)
//			{
//				delete [] cur_first_mapper_chrs->at(i_mapper_chr);
//				delete [] cur_first_mapper_read_ids->at(i_mapper_chr);
//			} // i_mapper_chr loop.
//			cur_first_mapper_chrs->clear();
//			cur_first_mapper_read_ids->clear();
//			
//			cur_last_mapper_indices->clear();
//			cur_last_mapper_strands->clear();
//			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_last_mapper_chrs->size(); i_mapper_chr++)
//			{
//				delete [] cur_last_mapper_chrs->at(i_mapper_chr);
//				delete [] cur_last_mapper_read_ids->at(i_mapper_chr);
//			} // i_mapper_chr loop.
//			cur_last_mapper_chrs->clear();
//			cur_last_mapper_read_ids->clear();
//			
//			// Add the currently read id.
//			if(cur_pair_index == 0)
//			{
//				cur_first_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
//				cur_first_mapper_indices->push_back(cur_start);
//				cur_first_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
//				cur_first_mapper_strands->push_back(cur_strand);
//			}
//			else
//			{
//				cur_last_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
//				cur_last_mapper_indices->push_back(cur_start);
//				cur_last_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
//				cur_last_mapper_strands->push_back(cur_strand);
//			}
//
//			// Copy the current read id.
//			t_string::copy(prev_read_id, cur_read_id);
//		} // prev_read id and cur_read_id comparison.
//
//		delete [] cur_read_line;
//	} // line reading loop.
//	fclose(f_sorted_pe_reads);
//
//	// Dump the connections between the bins.
//	FILE* f_pe_connections = open_f("pe_connections.list", "w");
//	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
//	{
//		for(int i_bin = 0; i_bin < (int)bin_regions_per_chr->at(i_chr)->size(); i_bin++)
//		{
//			// Dump the connections for the current region.
//			//vector<t_annot_region*>* cur_bin_connections = (vector<t_annot_region*>*)(bin_regions_per_chr->at(i_chr)->at(i_bin)->data);
//			vector<t_connection_info_entry*>* cur_bin_connections = (vector<t_connection_info_entry*>*)(bin_regions_per_chr->at(i_chr)->at(i_bin)->data);
//
//			// Process all the connections.
//			for(int i_conn = 0; i_conn < (int)cur_bin_connections->size(); i_conn++)
//			{
//				fprintf(f_pe_connections, "%s\t%d\t%d\t", bin_regions_per_chr->at(i_chr)->at(i_bin)->chrom, 
//					bin_regions_per_chr->at(i_chr)->at(i_bin)->start, bin_regions_per_chr->at(i_chr)->at(i_bin)->end);
//
//				fprintf(f_pe_connections, "%s\t%d\t%d\t%d\t", 
//					cur_bin_connections->at(i_conn)->connecting_region->chrom, 
//					cur_bin_connections->at(i_conn)->connecting_region->start, 
//					cur_bin_connections->at(i_conn)->connecting_region->end,
//					cur_bin_connections->at(i_conn)->read_ids->size());
//
//				for(int i_read = 0; i_read < cur_bin_connections->at(i_conn)->read_ids->size(); i_read++)
//				{
//					fprintf(f_pe_connections, "%s(%c, %c) \t", cur_bin_connections->at(i_conn)->read_ids->at(i_read), 
//						cur_bin_connections->at(i_conn)->first_read_strands->at(i_read), 
//						cur_bin_connections->at(i_conn)->last_read_strands->at(i_read));
//				} // i_read loop.
//
//				fprintf(f_pe_connections, "\n");
//
//			} // i_conn loop.
//		} // i_bin loop.
//	} // i_chr loop.
//
//	fclose(f_pe_connections);
//}

vector<int>* get_chromosome_lengths_per_mapped_reads(char* mapped_reads_dir)
{
	vector<int>* chr_lengths = new vector<int>();

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		int l_cur_chr = 0;
		char cur_line[1000];
		char cur_mapped_reads_fp[1000];
		sprintf(cur_mapped_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, chr_ids->at(i_chr));
		FILE* f_mapped_reads = get_processed_reads_ptr_wrapper(cur_mapped_reads_fp);

		while(1)
		{
			if(fgets(cur_line, 1000, f_mapped_reads) == NULL)
			{
				break;
			}

			int cur_pos = 0;
			sscanf(cur_line, "%*s %*s %d", &cur_pos);

			if(cur_pos > l_cur_chr)
			{
				l_cur_chr = cur_pos + 1000;
			}
		} // file reading loop.
		//fclose(f_mapped_reads);
		close_processed_reads_ptr_wrapper(f_mapped_reads, cur_mapped_reads_fp);

		chr_lengths->push_back(l_cur_chr);
	} // i_chr loop.

	return(chr_lengths);
}

void generate_discordant_PE_read_connections_per_sorted_PE_reads_file(char* chr_ids_lengths_fp, char* sorted_pe_reads_fp, 
	char* mapability_dir,
	int l_bin, 
	int l_step,
	int min_l_same_chr_separation, 
	int max_concordant_mapping_separation,
	char valid_first_mapper_strand,
	char valid_last_mapper_strand)
{
	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();

	load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_fp, chr_ids, chr_lengths);

	// Allocate the bins in each chromosome.
	vector<vector<t_annot_region*>*>* bin_regions_per_chr = new vector<vector<t_annot_region*>*>();
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_chr_mapp_prof_fp[1000];
		sprintf(cur_chr_mapp_prof_fp, "%s/%s.bin", mapability_dir, chr_ids->at(i_chr));
		int l_mapp_prof = 0;
		double* cur_chr_mapp_prof = NULL;
		if(check_file(cur_chr_mapp_prof_fp))
		{
			cur_chr_mapp_prof = load_per_nucleotide_binary_profile(cur_chr_mapp_prof_fp, l_mapp_prof);
		}
		else
		{
			fprintf(stderr, "Could not find mappability map %s\n", cur_chr_mapp_prof_fp);
		}

		fprintf(stderr, "Allocating bins in %s\n", chr_ids->at(i_chr));
		vector<t_annot_region*>* bin_regions_per_cur_chr = new vector<t_annot_region*>();

		int cur_bin_start = 1;
		//while(cur_bin_start + l_bin <= L_CHROM)
		while(cur_bin_start <= chr_lengths->at(i_chr))
		{
			double cur_reg_mapp = 5.0;
			if(cur_chr_mapp_prof != NULL)
			{
				cur_reg_mapp = 0;
				for(int i = cur_bin_start; i < cur_bin_start + l_bin; i++)
				{
					if(i < l_mapp_prof)
					{
						cur_reg_mapp += cur_chr_mapp_prof[i];
					}
				} // i loop.

				cur_reg_mapp /= l_bin;
			}

			t_annot_region* cur_bin_region = get_empty_region();
			cur_bin_region->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
			cur_bin_region->start = cur_bin_start;
			cur_bin_region->end = cur_bin_start + l_bin;
			t_per_bin_connection_info* cur_reg_conn_info = new t_per_bin_connection_info();
			cur_reg_conn_info->avg_mappability = cur_reg_mapp;
			cur_reg_conn_info->connecting_bins = new vector<t_connection_info_entry*>();
			cur_bin_region->data = cur_reg_conn_info;

			bin_regions_per_cur_chr->push_back(cur_bin_region);

			// Update bin start.
			cur_bin_start += l_step;
		} // i_bin loop.

		fprintf(stderr, "%ld regions (%d)\n", bin_regions_per_cur_chr->size(), cur_bin_start);
		bin_regions_per_chr->push_back(bin_regions_per_cur_chr);

		// Free mapp profile memory.
		if(cur_chr_mapp_prof != NULL)
		{
			delete [] cur_chr_mapp_prof;
		}
	} // i_chr loop.

	FILE* f_sorted_pe_reads = open_f(sorted_pe_reads_fp, "r");
	char prev_read_id[1000];
	prev_read_id[0] = 0;

	// Following is the mapping information for first and last mappers for the current read id.
	vector<char*>* cur_first_mapper_read_ids = new vector<char*>();
	vector<int>* cur_first_mapper_indices = new vector<int>();
	vector<char*>* cur_first_mapper_chrs = new vector<char*>();
	vector<char>* cur_first_mapper_strands = new vector<char>();
	vector<char*>* cur_last_mapper_read_ids = new vector<char*>();
	vector<int>* cur_last_mapper_indices = new vector<int>();
	vector<char*>* cur_last_mapper_chrs = new vector<char*>();
	vector<char>* cur_last_mapper_strands = new vector<char>();

	// Get the current list of connections: Get the consecutive lines with same read id.
	// MARILYN_0005:1:2:1272:993#0 76M R 102495555 0   5
	int n_lines_processed = 0;
	int n_concordant = 0;
	int n_discordant = 0;
	while(1)
	{
		if(n_lines_processed % 1000000 == 0)
		{
			fprintf(stderr, "Processed %d (%d, %d) read.           \r", n_lines_processed, n_concordant, n_discordant);
		}

		// Read the current read line.
		char* cur_read_line = getline(f_sorted_pe_reads);
		if(cur_read_line == NULL)
		{
			break;
		}

		n_lines_processed++;

		char cur_read_id[1000];
		char cur_mapping_str[1000];
		char cur_strand;
		int cur_start;
		int cur_pair_index;
		char cur_chr_id[100];

		if(sscanf(cur_read_line, "%s %s %c %d %d %s", 
			cur_read_id, cur_mapping_str, &cur_strand, &cur_start, &cur_pair_index, cur_chr_id) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", cur_read_line);
			exit(0);
		}

		// If this chromosome is not in the list, skip it.
		int i_chr = t_string::get_i_str(chr_ids, cur_chr_id);
		if(i_chr == (int)(chr_ids->size()))
		{
			delete [] cur_read_line;
			continue;
		}

		// Does this read have the same id as previous ones?
		if(t_string::compare_strings(prev_read_id, cur_read_id))
		{
			if(cur_pair_index == 0)
			{
				cur_first_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_first_mapper_indices->push_back(cur_start);
				cur_first_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_first_mapper_strands->push_back(cur_strand);
			}
			else
			{
				cur_last_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_last_mapper_indices->push_back(cur_start);
				cur_last_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_last_mapper_strands->push_back(cur_strand);
			}
		}
		else
		{
			bool found_concordant_mapping = false;
			for(int i_first = 0; 
				!found_concordant_mapping && i_first < (int)cur_first_mapper_indices->size(); 
				i_first++)
			{
				for(int i_last = 0; 
					!found_concordant_mapping && i_last < (int)cur_last_mapper_indices->size(); 
					i_last++)
				{
					// Update the current pair.
					char first_mapper_strand = cur_first_mapper_strands->at(i_first);
					char last_mapper_strand = cur_last_mapper_strands->at(i_last);
					char* first_mapper_chr = cur_first_mapper_chrs->at(i_first);
					char* last_mapper_chr = cur_last_mapper_chrs->at(i_last);

					// Is this the concordant mapping for this current read?
					if(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
						abs(cur_first_mapper_indices->at(i_first) - cur_last_mapper_indices->at(i_last)) < max_concordant_mapping_separation &&
						((first_mapper_strand == valid_first_mapper_strand && last_mapper_strand == valid_last_mapper_strand) ||
						(first_mapper_strand == valid_last_mapper_strand && last_mapper_strand == valid_first_mapper_strand)))
					{
						found_concordant_mapping = true;
					}
				} // i_last loop.
			} // i_first loop.

			/*********************************************************************************************************
			TODO: Classify the discordant reads with respect to the type of discordancy that identifies type of SV:
			.. 1.Distance
			.. 2.Strand
			.. 3.Mapping/Unapping
			*********************************************************************************************************/

			// Process the curent list of entries if there was no concordant mapping.
			if(!found_concordant_mapping)
			{
				n_discordant++;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
                if(cur_first_mapper_indices->size() > 0 && cur_last_mapper_indices->size() > 0)
                {
					fprintf(stderr, "First reads:\n");
					for(int i_first = 0; 
						i_first < (int)cur_first_mapper_indices->size(); 
						i_first++)
					{
						fprintf(stderr, "%s: %s: %d (%c)\n", cur_first_mapper_read_ids->at(i_first), 
							cur_first_mapper_chrs->at(i_first), 
							cur_first_mapper_indices->at(i_first), 
							cur_first_mapper_strands->at(i_first));
					} // i_first loop.

					fprintf(stderr, "Last reads:\n");
					for(int i_last = 0; 
						i_last < (int)cur_last_mapper_indices->size(); 
						i_last++)
					{
						fprintf(stderr, "%s: %s: %d (%c)\n", cur_last_mapper_read_ids->at(i_last), 
							cur_last_mapper_chrs->at(i_last),
							cur_last_mapper_indices->at(i_last), 
							cur_last_mapper_strands->at(i_last));
					} // i_last loop.
					getc(stdin);
				}
} // check for dumping reads.

				/****************************************************************************************
				TODO: Rather than reporting all the discorant pairs, use the pair that has the highest
				mapping qualities.
				****************************************************************************************/
				for(int i_first = 0; 
					i_first < (int)cur_first_mapper_indices->size(); 
					i_first++)
				{
					for(int i_last = 0; 
						i_last < (int)cur_last_mapper_indices->size(); 
						i_last++)
					{
						// Update the current pair.
						char first_mapper_strand = cur_first_mapper_strands->at(i_first);
						char last_mapper_strand = cur_last_mapper_strands->at(i_last);

						int first_mapper_index = cur_first_mapper_indices->at(i_first);
						int last_mapper_index = cur_last_mapper_indices->at(i_last);

						//char* first_mapper_chr = cur_first_mapper_chrs->at(i_first);
						//char* last_mapper_chr = cur_last_mapper_chrs->at(i_last);

						// Since we know that the pairs are discordant, we do not need any more filters.
						//if(first_mapper_strand != last_mapper_strand &&			
						//	(!t_string::compare_strings(first_mapper_chr, last_mapper_chr) || 
						//	(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
						//	abs(cur_first_mapper_indices->at(i_first) - cur_last_mapper_indices->at(i_last)) > min_l_same_chr_separation)))
						//if(!t_string::compare_strings(first_mapper_chr, last_mapper_chr))
						{
							// Get the chr indices for the mapping target bins.
							int i_first_chr_id = t_string::get_i_str(chr_ids, cur_first_mapper_chrs->at(i_first));
							int i_last_chr_id = t_string::get_i_str(chr_ids, cur_last_mapper_chrs->at(i_last));

							if(i_first_chr_id == (int)chr_ids->size() || 
								i_last_chr_id == (int)chr_ids->size())
							{
								fprintf(stderr, "Could not locate the chromosome id.\n");
								exit(0);
							}

							t_annot_region* cur_first_bin = bin_regions_per_chr->at(i_first_chr_id)->at(cur_first_mapper_indices->at(i_first) / l_bin);
							t_annot_region* cur_last_bin = bin_regions_per_chr->at(i_last_chr_id)->at(cur_last_mapper_indices->at(i_last) / l_bin);

							vector<t_connection_info_entry*>* cur_first_bin_bins_list = ((t_per_bin_connection_info*)(vector<t_connection_info_entry*>*)(cur_first_bin->data))->connecting_bins;
							vector<t_connection_info_entry*>* cur_last_bin_bins_list = ((t_per_bin_connection_info*)(vector<t_connection_info_entry*>*)(cur_last_bin->data))->connecting_bins;

							// Add the current pair to the list of connections.
							bool found_connection1 = false;
							for(int i_conn1 = 0; 
								i_conn1 < (int)cur_first_bin_bins_list->size() &&
								!found_connection1; 
								i_conn1++)
							{
								if(t_string::compare_strings(cur_last_bin->chrom, cur_first_bin_bins_list->at(i_conn1)->connecting_region->chrom) &&
									cur_last_bin->start == cur_first_bin_bins_list->at(i_conn1)->connecting_region->start &&
									cur_last_bin->end == cur_first_bin_bins_list->at(i_conn1)->connecting_region->end)
								{
									cur_first_bin_bins_list->at(i_conn1)->read_ids->push_back(t_string::copy_me_str(prev_read_id));
									cur_first_bin_bins_list->at(i_conn1)->first_read_strands->push_back(first_mapper_strand);
									cur_first_bin_bins_list->at(i_conn1)->first_read_starts->push_back(first_mapper_index);
									cur_first_bin_bins_list->at(i_conn1)->last_read_strands->push_back(last_mapper_strand);
									cur_first_bin_bins_list->at(i_conn1)->last_read_starts->push_back(last_mapper_index);

									cur_first_bin_bins_list->at(i_conn1)->first_read_i_chrs->push_back(i_first_chr_id);
									cur_first_bin_bins_list->at(i_conn1)->last_read_i_chrs->push_back(i_last_chr_id);
									
									found_connection1 = true;
								}
							} // i_conn1 loop.

							if(!found_connection1)
							{
								t_connection_info_entry* new_first_connection_entry = new t_connection_info_entry();
								new_first_connection_entry->connecting_region = cur_last_bin;
								new_first_connection_entry->read_ids = new vector<char*>();
								new_first_connection_entry->read_ids->push_back(t_string::copy_me_str(prev_read_id));
								new_first_connection_entry->first_read_strands = new vector<char>();
								new_first_connection_entry->first_read_strands->push_back(first_mapper_strand);
								new_first_connection_entry->last_read_strands = new vector<char>();
								new_first_connection_entry->last_read_strands->push_back(last_mapper_strand);
								new_first_connection_entry->first_read_starts = new vector<int>();
								new_first_connection_entry->first_read_starts->push_back(first_mapper_index);
								new_first_connection_entry->last_read_starts = new vector<int>();
								new_first_connection_entry->last_read_starts->push_back(last_mapper_index);
								new_first_connection_entry->first_read_i_chrs = new vector<int>();
								new_first_connection_entry->first_read_i_chrs->push_back(i_first_chr_id);
								new_first_connection_entry->last_read_i_chrs = new vector<int>();
								new_first_connection_entry->last_read_i_chrs->push_back(i_last_chr_id);

								cur_first_bin_bins_list->push_back(new_first_connection_entry);
							}

							// Look for the first bin in the connection list of last bin.
							// If the bins are the same bin, do not run following as it is redundant.
							bool found_connection2 = false;
							if(cur_first_bin != cur_last_bin)
							{
								for(int i_conn1 = 0; 
									i_conn1 < (int)cur_last_bin_bins_list->size() &&
									!found_connection2; 
									i_conn1++)
								{
									if(t_string::compare_strings(cur_first_bin->chrom, cur_last_bin_bins_list->at(i_conn1)->connecting_region->chrom) &&
										cur_first_bin->start ==  cur_last_bin_bins_list->at(i_conn1)->connecting_region->start &&
										cur_first_bin->end == cur_last_bin_bins_list->at(i_conn1)->connecting_region->end)
									{
										cur_last_bin_bins_list->at(i_conn1)->read_ids->push_back(t_string::copy_me_str(prev_read_id));
										cur_last_bin_bins_list->at(i_conn1)->first_read_strands->push_back(first_mapper_strand);
										cur_last_bin_bins_list->at(i_conn1)->last_read_strands->push_back(last_mapper_strand);
										cur_last_bin_bins_list->at(i_conn1)->first_read_starts->push_back(first_mapper_index);
										cur_last_bin_bins_list->at(i_conn1)->last_read_starts->push_back(last_mapper_index);

										cur_last_bin_bins_list->at(i_conn1)->first_read_i_chrs->push_back(i_first_chr_id);
										cur_last_bin_bins_list->at(i_conn1)->last_read_i_chrs->push_back(i_last_chr_id);
										found_connection2 = true;
									}
								} // i_conn1 loop.

								if(!found_connection2)
								{
									t_connection_info_entry* new_last_connection_entry = new t_connection_info_entry();							
									new_last_connection_entry->connecting_region = cur_first_bin;
									new_last_connection_entry->read_ids = new vector<char*>();
									new_last_connection_entry->read_ids->push_back(t_string::copy_me_str(prev_read_id));
									new_last_connection_entry->first_read_strands = new vector<char>();
									new_last_connection_entry->first_read_strands->push_back(first_mapper_strand);
									new_last_connection_entry->last_read_strands = new vector<char>();
									new_last_connection_entry->last_read_strands->push_back(last_mapper_strand);
									new_last_connection_entry->first_read_starts = new vector<int>();
									new_last_connection_entry->first_read_starts->push_back(first_mapper_index);
									new_last_connection_entry->last_read_starts = new vector<int>();
									new_last_connection_entry->last_read_starts->push_back(last_mapper_index);

									new_last_connection_entry->first_read_i_chrs = new vector<int>();
									new_last_connection_entry->first_read_i_chrs->push_back(i_first_chr_id);
									new_last_connection_entry->last_read_i_chrs = new vector<int>();
									new_last_connection_entry->last_read_i_chrs->push_back(i_last_chr_id);

									cur_last_bin_bins_list->push_back(new_last_connection_entry);
								}
							} // self check for cutting this loop.

							if(cur_first_bin != cur_last_bin &&
								found_connection2 != found_connection1)
							{
								fprintf(stderr, "The connections are not consistent: %s\t%d\t%d, %s\t%d\t%d\nLeft mapping: %d(%c), Right mapping: %d(%c)", 
									cur_first_bin->chrom, cur_first_bin->start, cur_first_bin->end, 
									cur_last_bin->chrom, cur_last_bin->start, cur_last_bin->end,
									cur_first_mapper_indices->at(i_first), first_mapper_strand, 
									cur_last_mapper_indices->at(i_last), last_mapper_strand);
								exit(0);
							}
						} // extra discordance checks.
					} // i_last loop.
				} // i_first loop.
			} // concordant read check.
			else
			{
				n_concordant++;
			}

			// Process the current entry, clear the list, then add the entry for the current read.
			cur_first_mapper_indices->clear();
			cur_first_mapper_strands->clear();
			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_first_mapper_chrs->size(); i_mapper_chr++)
			{
				delete [] cur_first_mapper_chrs->at(i_mapper_chr);
				delete [] cur_first_mapper_read_ids->at(i_mapper_chr);
			} // i_mapper_chr loop.
			cur_first_mapper_chrs->clear();
			cur_first_mapper_read_ids->clear();
			
			cur_last_mapper_indices->clear();
			cur_last_mapper_strands->clear();
			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_last_mapper_chrs->size(); i_mapper_chr++)
			{
				delete [] cur_last_mapper_chrs->at(i_mapper_chr);
				delete [] cur_last_mapper_read_ids->at(i_mapper_chr);
			} // i_mapper_chr loop.
			cur_last_mapper_chrs->clear();
			cur_last_mapper_read_ids->clear();
			
			// Add the currently read id.
			if(cur_pair_index == 0)
			{
				cur_first_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_first_mapper_indices->push_back(cur_start);
				cur_first_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_first_mapper_strands->push_back(cur_strand);
			}
			else
			{
				cur_last_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_last_mapper_indices->push_back(cur_start);
				cur_last_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_last_mapper_strands->push_back(cur_strand);
			}

			// Copy the current read id.
			t_string::copy(prev_read_id, cur_read_id);
		} // prev_read id and cur_read_id comparison.

		delete [] cur_read_line;
	} // line reading loop.
	fclose(f_sorted_pe_reads);

	// Dump the connections between the bins.
	FILE* f_pe_connections = open_f("pe_connections.list", "w");
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_chr_mapp_prof_fp[1000];
		sprintf(cur_chr_mapp_prof_fp, "%s/%s.bin", mapability_dir, chr_ids->at(i_chr));
		int l_mapp_prof = 0;
		double* cur_chr_mapp_prof = NULL;
		if(check_file(cur_chr_mapp_prof_fp))
		{
			cur_chr_mapp_prof = load_per_nucleotide_binary_profile(cur_chr_mapp_prof_fp, l_mapp_prof);
		}
		else
		{
			fprintf(stderr, "Could not find mappability map %s\n", cur_chr_mapp_prof_fp);
		}

		for(int i_bin = 0; i_bin < (int)bin_regions_per_chr->at(i_chr)->size(); i_bin++)
		{
			// Dump the connections for the current region.
			vector<t_connection_info_entry*>* cur_bin_connections = ((t_per_bin_connection_info*)(vector<t_connection_info_entry*>*)(bin_regions_per_chr->at(i_chr)->at(i_bin)->data))->connecting_bins;
			double avg_mapp = ((t_per_bin_connection_info*)(vector<t_connection_info_entry*>*)(bin_regions_per_chr->at(i_chr)->at(i_bin)->data))->avg_mappability;

			// Process all the connections.
			for(int i_conn = 0; i_conn < (int)cur_bin_connections->size(); i_conn++)
			{
				double avg_read_start_mapp = 12345;
				if(cur_chr_mapp_prof != NULL)
				{
					avg_read_start_mapp = 0;
					for(int i_read = 0; i_read < (int)cur_bin_connections->at(i_conn)->read_ids->size(); i_read++)
					{
						// One of the side's chromosome has to match the current chromosome.
						if(i_chr == cur_bin_connections->at(i_conn)->first_read_i_chrs->at(i_read))
						{
							avg_read_start_mapp += cur_chr_mapp_prof[cur_bin_connections->at(i_conn)->first_read_starts->at(i_read)];
						}
						else if(i_chr == cur_bin_connections->at(i_conn)->last_read_i_chrs->at(i_read))
						{
							avg_read_start_mapp += cur_chr_mapp_prof[cur_bin_connections->at(i_conn)->last_read_starts->at(i_read)];
						}
						else
						{
							fprintf(stderr, "Sanity check failed: %s(%d)\n", __FILE__, __LINE__);
							exit(0);
						}
					} // i_read loop.

					// Take the average of the read start mapabilities.
					if(cur_bin_connections->at(i_conn)->read_ids->size() > 0)
					{
						avg_read_start_mapp /= cur_bin_connections->at(i_conn)->read_ids->size();
					}
					else
					{
						fprintf(stderr, "Sanity check failed: %s(%d)\n", __FILE__, __LINE__);
						exit(0);
					}
				} // mappability map check.

				// Dump the info.
				fprintf(f_pe_connections, "%s\t%d\t%d\t%lf\t%lf\t", bin_regions_per_chr->at(i_chr)->at(i_bin)->chrom, 
					bin_regions_per_chr->at(i_chr)->at(i_bin)->start, bin_regions_per_chr->at(i_chr)->at(i_bin)->end,
					avg_mapp, avg_read_start_mapp);

				fprintf(f_pe_connections, "%s\t%d\t%d\t%d\t", 
					cur_bin_connections->at(i_conn)->connecting_region->chrom, 
					cur_bin_connections->at(i_conn)->connecting_region->start, 
					cur_bin_connections->at(i_conn)->connecting_region->end,
					(int)cur_bin_connections->at(i_conn)->read_ids->size());

				for(int i_read = 0; i_read < (int)cur_bin_connections->at(i_conn)->read_ids->size(); i_read++)
				{
					fprintf(f_pe_connections, "%s(%c, %c) \t", cur_bin_connections->at(i_conn)->read_ids->at(i_read), 
						cur_bin_connections->at(i_conn)->first_read_strands->at(i_read), 
						cur_bin_connections->at(i_conn)->last_read_strands->at(i_read));
				} // i_read loop.

				fprintf(f_pe_connections, "\n");
			} // i_conn loop.
		} // i_bin loop.

		// Free mapp profile memory.
		if(cur_chr_mapp_prof != NULL)
		{
			delete [] cur_chr_mapp_prof;
		}
	} // i_chr loop.

	fclose(f_pe_connections);
}

void generate_discordant_read_connections_per_id_sorted_SE_mapped_reads_file(char* chr_ids_fp, char* sorted_se_mapped_reads_fp, 
	int l_bin, int min_l_same_chr_separation, int max_concordant_mapping_separation,
	char valid_first_mapper_strand,
	char valid_last_mapper_strand)
{
	fprintf(stderr, "generate_discordant_read_connections_per_id_sorted_SE_mapped_reads_file is currently conceptually defunct and should not be used. Exiting.\n");
	exit(0);

//	// Allocate the bins for counting the connections.
//	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
//	int n_bins_per_chr = (L_CHROM / l_bin) + 10;
//
//	// Allocate the bins in each chromosome.
//	vector<vector<t_annot_region*>*>* bin_regions_per_chr = new vector<vector<t_annot_region*>*>();
//	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
//	{		
//		fprintf(stderr, "Allocating bins in %s\n", chr_ids->at(i_chr));
//		vector<t_annot_region*>* bin_regions_per_cur_chr = new vector<t_annot_region*>();
//
//		int cur_bin_start = 1;
//		for(int i_bin = 0; i_bin < n_bins_per_chr; i_bin++)
//		{
//			t_annot_region* cur_bin_region = get_empty_region();
//			cur_bin_region->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
//			cur_bin_region->start = cur_bin_start;
//			cur_bin_region->end = cur_bin_start + l_bin;
//			vector<t_connection_info_entry*>* cur_bin_connecting_bins = new vector<t_connection_info_entry*>();
//			cur_bin_region->data = cur_bin_connecting_bins;
//
//			bin_regions_per_cur_chr->push_back(cur_bin_region);
//
//			cur_bin_start += l_bin;
//		} // i_bin loop.
//
//		fprintf(stderr, "%ld regions (%d)\n", bin_regions_per_cur_chr->size(), cur_bin_start);
//		bin_regions_per_chr->push_back(bin_regions_per_cur_chr);
//	} // i_chr loop.
//
//	FILE* f_sorted_pe_reads = open_f(sorted_se_mapped_reads_fp, "r");
//	char prev_read_id[1000];
//	prev_read_id[0] = 0;
//
//	vector<char*>* cur_read_ids = new vector<char*>();
//	vector<int>* cur_indices = new vector<int>();
//	vector<char*>* cur_chrs = new vector<char*>();
//	vector<char>* cur_strands = new vector<char>();
//
//	// Get the current list of connections: Get the consecutive lines with same read id.
//	// MARILYN_0005:1:2:1272:993#0 76M R 102495555 0   5
//	int n_lines_processed = 0;
//	while(1)
//	{
//		if(n_lines_processed % 1000000 == 0)
//		{
//			fprintf(stderr, "Processed %d read.           \r", n_lines_processed);
//		}
//
//		// Read the current read line.
//		char* cur_read_line = getline(f_sorted_pe_reads);
//		if(cur_read_line == NULL)
//		{
//			break;
//		}
//
//		n_lines_processed++;
//
//		char cur_read_id[1000];
//		char cur_mapping_str[1000];
//		char cur_strand;
//		int cur_start;
//		int cur_pair_index;
//		char cur_chr_id[100];
//
//		if(sscanf(cur_read_line, "%s %s %c %d %d %s", 
//			cur_read_id, cur_mapping_str, &cur_strand, &cur_start, &cur_pair_index, cur_chr_id) != 6)
//		{
//			fprintf(stderr, "Could not parse: %s\n", cur_read_line);
//			exit(0);
//		}
//
//		// Does this read have the same id as previous ones?
//		if(t_string::compare_strings(prev_read_id, cur_read_id))
//		{
//			cur_read_ids->push_back(t_string::copy_me_str(cur_read_id));
//			cur_indices->push_back(cur_start);
//			cur_chrs->push_back(t_string::copy_me_str(cur_chr_id));
//			cur_strands->push_back(cur_strand);
//		}
//		else
//		{
//			bool found_concordant_mapping = false;
//			for(int i_first = 0; 
//				!found_concordant_mapping && i_first < (int)cur_indices->size(); 
//				i_first++)
//			{
//				for(int i_last = i_first+1;
//					!found_concordant_mapping && i_last < (int)cur_indices->size(); 
//					i_last++)
//				{
//					// Update the current pair.
//					char first_mapper_strand = cur_strands->at(i_first);
//					char last_mapper_strand = cur_strands->at(i_last);
//					char* first_mapper_chr = cur_chrs->at(i_first);
//					char* last_mapper_chr = cur_chrs->at(i_last);
//
//					// Is this the concordant mapping for this current read?
//					if(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
//						abs(cur_indices->at(i_first) - cur_indices->at(i_last)) < max_concordant_mapping_separation &&
//						((first_mapper_strand == valid_first_mapper_strand && last_mapper_strand == valid_last_mapper_strand) ||
//						(first_mapper_strand == valid_last_mapper_strand && last_mapper_strand == valid_first_mapper_strand)))
//					{
//						found_concordant_mapping = true;
//					}
//				} // i_last loop.
//			} // i_first loop.
//
//			// Process the curent list of entries if there was no concordant mapping.
//			if(!found_concordant_mapping)
//			{
//if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
//{
//                if(cur_indices->size() > 0 && cur_indices->size() > 0)
//                {
//					fprintf(stderr, "First reads:\n");
//					for(int i_first = 0; 
//						i_first < (int)cur_indices->size(); 
//						i_first++)
//					{
//						fprintf(stderr, "%s: %s: %d (%c)\n", cur_indices->at(i_first), 
//							cur_chrs->at(i_first), 
//							cur_indices->at(i_first), 
//							cur_strands->at(i_first));
//					} // i_first loop.
//
//					fprintf(stderr, "Last reads:\n");
//					for(int i_last = 0; 
//						i_last < (int)cur_indices->size(); 
//						i_last++)
//					{
//						fprintf(stderr, "%s: %s: %d (%c)\n", cur_read_ids->at(i_last), 
//							cur_chrs->at(i_last),
//							cur_indices->at(i_last), 
//							cur_strands->at(i_last));
//					} // i_last loop.
//					getc(stdin);
//				}
//} // check for dumping reads.
//
//				for(int i_first = 0; 
//					i_first < (int)cur_indices->size(); 
//					i_first++)
//				{
//					for(int i_last = i_first+1; 
//						i_last < (int)cur_indices->size(); 
//						i_last++)
//					{
//						// Update the current pair.
//						char first_mapper_strand = cur_strands->at(i_first);
//						char last_mapper_strand = cur_strands->at(i_last);
//
//						char* first_mapper_chr = cur_chrs->at(i_first);
//						char* last_mapper_chr = cur_chrs->at(i_last);
//
//						// Since we know that the pairs are discordant, we do not need any more filters.
//						//if(first_mapper_strand != last_mapper_strand &&			
//						//	(!t_string::compare_strings(first_mapper_chr, last_mapper_chr) || 
//						//	(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
//						//	abs(cur_first_mapper_indices->at(i_first) - cur_last_mapper_indices->at(i_last)) > min_l_same_chr_separation)))
//						//if(!t_string::compare_strings(first_mapper_chr, last_mapper_chr))
//						{
//							// Get the chr indices for the mapping target bins.
//							int i_first_chr_id = t_string::get_i_str(chr_ids, cur_chrs->at(i_first));
//							int i_last_chr_id = t_string::get_i_str(chr_ids, cur_chrs->at(i_last));
//
//							t_annot_region* cur_first_bin = bin_regions_per_chr->at(i_first_chr_id)->at(cur_indices->at(i_first) / l_bin);
//							t_annot_region* cur_last_bin = bin_regions_per_chr->at(i_last_chr_id)->at(cur_indices->at(i_last) / l_bin);
//
//							vector<t_connection_info_entry*>* cur_first_bin_bins_list = (vector<t_connection_info_entry*>*)(cur_first_bin->data);
//							vector<t_connection_info_entry*>* cur_last_bin_bins_list = (vector<t_connection_info_entry*>*)(cur_last_bin->data);
//
//							// Add the regions to the connection lists, if they do not exist, add them, otherwise include them.
//							bool found_connection1 = false;
//							for(int i_conn1 = 0; 
//								i_conn1 < cur_first_bin_bins_list->size() &&
//								!found_connection1; 
//								i_conn1++)
//							{
//								if(t_string::compare_strings(cur_last_bin->chrom, cur_first_bin_bins_list->at(i_conn1)->connecting_region->chrom) &&
//									cur_last_bin->start == cur_first_bin_bins_list->at(i_conn1)->connecting_region->start &&
//									cur_last_bin->end == cur_first_bin_bins_list->at(i_conn1)->connecting_region->end)
//								{
//									cur_first_bin_bins_list->at(i_conn1)->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//									cur_first_bin_bins_list->at(i_conn1)->first_read_strands->push_back(first_mapper_strand);
//									cur_first_bin_bins_list->at(i_conn1)->last_read_strands->push_back(last_mapper_strand);
//									found_connection1 = true;
//								}
//							} // i_conn1 loop.
//
//							if(!found_connection1)
//							{
//								t_connection_info_entry* new_first_connection_entry = new t_connection_info_entry();
//								new_first_connection_entry->connecting_region = cur_last_bin;
//								new_first_connection_entry->read_ids = new vector<char*>();
//								new_first_connection_entry->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//								new_first_connection_entry->first_read_strands = new vector<char>();
//								new_first_connection_entry->first_read_strands->push_back(first_mapper_strand);
//								new_first_connection_entry->last_read_strands = new vector<char>();
//								new_first_connection_entry->last_read_strands->push_back(last_mapper_strand);
//								cur_first_bin_bins_list->push_back(new_first_connection_entry);
//							}
//
//							// Look for the first bin in the connection list of last bin.
//							// If the bins are the same bin, do not run following as it is redundant.
//							bool found_connection2 = false;
//							if(cur_first_bin != cur_last_bin)
//							{
//								for(int i_conn1 = 0; 
//									i_conn1 < cur_last_bin_bins_list->size() &&
//									!found_connection2; 
//									i_conn1++)
//								{
//									if(t_string::compare_strings(cur_first_bin->chrom, cur_last_bin_bins_list->at(i_conn1)->connecting_region->chrom) &&
//										cur_first_bin->start ==  cur_last_bin_bins_list->at(i_conn1)->connecting_region->start &&
//										cur_first_bin->end == cur_last_bin_bins_list->at(i_conn1)->connecting_region->end)
//									{
//										cur_last_bin_bins_list->at(i_conn1)->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//										cur_last_bin_bins_list->at(i_conn1)->first_read_strands->push_back(first_mapper_strand);
//										cur_last_bin_bins_list->at(i_conn1)->last_read_strands->push_back(last_mapper_strand);
//										found_connection2 = true;
//									}
//								} // i_conn1 loop.
//
//								if(!found_connection2)
//								{
//									t_connection_info_entry* new_last_connection_entry = new t_connection_info_entry();							
//									new_last_connection_entry->connecting_region = cur_first_bin;
//									new_last_connection_entry->read_ids = new vector<char*>();
//									new_last_connection_entry->read_ids->push_back(t_string::copy_me_str(prev_read_id));
//									new_last_connection_entry->first_read_strands = new vector<char>();
//									new_last_connection_entry->first_read_strands->push_back(first_mapper_strand);
//									new_last_connection_entry->last_read_strands = new vector<char>();
//									new_last_connection_entry->last_read_strands->push_back(last_mapper_strand);
//									cur_last_bin_bins_list->push_back(new_last_connection_entry);
//								}
//							} // self check for cutting this loop.
//
//							if(cur_first_bin != cur_last_bin &&
//								found_connection2 != found_connection1)
//							{
//								fprintf(stderr, "The connections are not consistent: %s\t%d\t%d, %s\t%d\t%d\nLeft mapping: %d(%c), Right mapping: %d(%c)", 
//									cur_first_bin->chrom, cur_first_bin->start, cur_first_bin->end, 
//									cur_last_bin->chrom, cur_last_bin->start, cur_last_bin->end,
//									cur_indices->at(i_first), first_mapper_strand, 
//									cur_indices->at(i_last), last_mapper_strand);
//								exit(0);
//							}
//						} // extra discordance checks.
//					} // i_last loop.
//				} // i_first loop.
//			} // concordant read check.
//
//			// Process the current entry, clear the list, then add the entry for the current read.
//			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_chrs->size(); i_mapper_chr++)
//			{
//				delete [] cur_chrs->at(i_mapper_chr);
//				delete [] cur_read_ids->at(i_mapper_chr);
//			} // i_mapper_chr loop.
//			cur_chrs->clear();
//			cur_read_ids->clear();			
//			cur_indices->clear();
//			cur_strands->clear();
//			
//			// Add the currently read id.
//			cur_read_ids->push_back(t_string::copy_me_str(cur_read_id));
//			cur_indices->push_back(cur_start);
//			cur_chrs->push_back(t_string::copy_me_str(cur_chr_id));
//			cur_strands->push_back(cur_strand);
//
//			// Copy the current read id.
//			t_string::copy(prev_read_id, cur_read_id);
//		} // prev_read id and cur_read_id comparison.
//
//		delete [] cur_read_line;
//	} // line reading loop.
//	fclose(f_sorted_pe_reads);
//
//	// Dump the connections between the bins.
//	FILE* f_pe_connections = open_f("pe_connections.list", "w");
//	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
//	{
//		for(int i_bin = 0; i_bin < (int)bin_regions_per_chr->at(i_chr)->size(); i_bin++)
//		{
//			// Dump the connections for the current region.
//			//vector<t_annot_region*>* cur_bin_connections = (vector<t_annot_region*>*)(bin_regions_per_chr->at(i_chr)->at(i_bin)->data);
//			vector<t_connection_info_entry*>* cur_bin_connections = (vector<t_connection_info_entry*>*)(bin_regions_per_chr->at(i_chr)->at(i_bin)->data);
//
//			// Process all the connections.
//			for(int i_conn = 0; i_conn < (int)cur_bin_connections->size(); i_conn++)
//			{
//				fprintf(f_pe_connections, "%s\t%d\t%d\t", bin_regions_per_chr->at(i_chr)->at(i_bin)->chrom, 
//					bin_regions_per_chr->at(i_chr)->at(i_bin)->start, bin_regions_per_chr->at(i_chr)->at(i_bin)->end);
//
//				fprintf(f_pe_connections, "%s\t%d\t%d\t%d\t", 
//					cur_bin_connections->at(i_conn)->connecting_region->chrom, 
//					cur_bin_connections->at(i_conn)->connecting_region->start, 
//					cur_bin_connections->at(i_conn)->connecting_region->end,
//					cur_bin_connections->at(i_conn)->read_ids->size());
//
//				for(int i_read = 0; i_read < cur_bin_connections->at(i_conn)->read_ids->size(); i_read++)
//				{
//					fprintf(f_pe_connections, "%s(%c, %c) \t", cur_bin_connections->at(i_conn)->read_ids->at(i_read), 
//						cur_bin_connections->at(i_conn)->first_read_strands->at(i_read), 
//						cur_bin_connections->at(i_conn)->last_read_strands->at(i_read));
//				} // i_read loop.
//
//				fprintf(f_pe_connections, "\n");
//			} // i_conn loop.
//		} // i_bin loop.
//	} // i_chr loop.
//
//	fclose(f_pe_connections);
} // generate_discordant_read_connections_per_id_sorted_SE_mapped_reads_file

void get_insert_stats_per_concordant_PE_reads(char* sorted_pe_reads_fp, 
	int l_read,
	int max_concordant_mapping_separation, 
	char valid_first_mapper_strand,
	char valid_last_mapper_strand)
{
	FILE* f_sorted_pe_reads = open_f(sorted_pe_reads_fp, "r");
	char prev_read_id[1000];
	prev_read_id[0] = 0;

	vector<char*>* cur_first_mapper_read_ids = new vector<char*>();
	vector<int>* cur_first_mapper_indices = new vector<int>();
	vector<char*>* cur_first_mapper_chrs = new vector<char*>();
	vector<char>* cur_first_mapper_strands = new vector<char>();
	vector<char*>* cur_last_mapper_read_ids = new vector<char*>();
	vector<int>* cur_last_mapper_indices = new vector<int>();
	vector<char*>* cur_last_mapper_chrs = new vector<char*>();
	vector<char>* cur_last_mapper_strands = new vector<char>();

	// Get the current list of connections: Get the consecutive lines with same read id.
	// MARILYN_0005:1:2:1272:993#0 76M R 102495555 0   5
	vector<int>* insert_l_per_mappings = new vector<int>();
	int n_lines_processed = 0;
	while(1)
	{
		if(n_lines_processed % 1000000 == 0)
		{
			fprintf(stderr, "Processed %d read.           \r", n_lines_processed);
		}

		// Read the current read line.
		char* cur_read_line = getline(f_sorted_pe_reads);
		if(cur_read_line == NULL)
		{
			break;
		}

		n_lines_processed++;

		char cur_read_id[1000];
		char cur_mapping_str[1000];
		char cur_strand;
		int cur_start;
		int cur_pair_index;
		char cur_chr_id[100];

		if(sscanf(cur_read_line, "%s %s %c %d %d %s", 
			cur_read_id, cur_mapping_str, &cur_strand, &cur_start, &cur_pair_index, cur_chr_id) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", cur_read_line);
			exit(0);
		}

		// Does this read have the same id as previous ones?
		if(t_string::compare_strings(prev_read_id, cur_read_id))
		{
			if(cur_pair_index == 0)
			{
				cur_first_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_first_mapper_indices->push_back(cur_start);
				cur_first_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_first_mapper_strands->push_back(cur_strand);
			}
			else
			{
				cur_last_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_last_mapper_indices->push_back(cur_start);
				cur_last_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_last_mapper_strands->push_back(cur_strand);
			}
		}
		else
		{
			for(int i_first = 0; 
				i_first < (int)cur_first_mapper_indices->size(); 
				i_first++)
			{
				for(int i_last = 0; 
					i_last < (int)cur_last_mapper_indices->size(); 
					i_last++)
				{
					// Update the current pair.
					char first_mapper_strand = cur_first_mapper_strands->at(i_first);
					char last_mapper_strand = cur_last_mapper_strands->at(i_last);
					char* first_mapper_chr = cur_first_mapper_chrs->at(i_first);
					char* last_mapper_chr = cur_last_mapper_chrs->at(i_last);

					// Is this the concordant mapping for this current read?
					if(t_string::compare_strings(first_mapper_chr, last_mapper_chr) && 
						abs(cur_first_mapper_indices->at(i_first) - cur_last_mapper_indices->at(i_last)) < max_concordant_mapping_separation &&
						((first_mapper_strand == valid_first_mapper_strand && last_mapper_strand == valid_last_mapper_strand) ||
						(first_mapper_strand == valid_last_mapper_strand && last_mapper_strand == valid_first_mapper_strand)))
					{
						int cur_insert_l = abs(cur_last_mapper_indices->at(i_last) - cur_first_mapper_indices->at(i_first)) + l_read;
						insert_l_per_mappings->push_back(cur_insert_l);
					}
				} // i_last loop.
			} // i_first loop.

			// Process the current entry, clear the list, then add the entry for the current read.
			cur_first_mapper_indices->clear();
			cur_first_mapper_strands->clear();
			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_first_mapper_chrs->size(); i_mapper_chr++)
			{
				delete [] cur_first_mapper_chrs->at(i_mapper_chr);
				delete [] cur_first_mapper_read_ids->at(i_mapper_chr);
			} // i_mapper_chr loop.
			cur_first_mapper_chrs->clear();
			cur_first_mapper_read_ids->clear();
			
			cur_last_mapper_indices->clear();
			cur_last_mapper_strands->clear();
			for(int i_mapper_chr = 0; i_mapper_chr < (int)cur_last_mapper_chrs->size(); i_mapper_chr++)
			{
				delete [] cur_last_mapper_chrs->at(i_mapper_chr);
				delete [] cur_last_mapper_read_ids->at(i_mapper_chr);
			} // i_mapper_chr loop.
			cur_last_mapper_chrs->clear();
			cur_last_mapper_read_ids->clear();
			
			// Add the currently read id.
			if(cur_pair_index == 0)
			{
				cur_first_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_first_mapper_indices->push_back(cur_start);
				cur_first_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_first_mapper_strands->push_back(cur_strand);
			}
			else
			{
				cur_last_mapper_read_ids->push_back(t_string::copy_me_str(cur_read_id));
				cur_last_mapper_indices->push_back(cur_start);
				cur_last_mapper_chrs->push_back(t_string::copy_me_str(cur_chr_id));
				cur_last_mapper_strands->push_back(cur_strand);
			}

			// Copy the current read id.
			t_string::copy(prev_read_id, cur_read_id);
		} // prev_read id and cur_read_id comparison.

		delete [] cur_read_line;
	} // line reading loop.
	fclose(f_sorted_pe_reads);

	// Get the statistics.
	sort(insert_l_per_mappings->begin(), insert_l_per_mappings->end());
	fprintf(stderr, "Median insert size is %d\n", insert_l_per_mappings->at(insert_l_per_mappings->size() / 2));
}

void get_differential_discordant_PE_pair_connections(char* sample_fp, char* control_fp, char* op_fp)
{
	// Load the sample file.
	vector<t_annot_region*>* sample_pe_regions = load_BED_with_line_information(sample_fp);
	if(sample_pe_regions == NULL)
	{
		fprintf(stderr, "Could not load %s\n", sample_fp);
		exit(0);
	}

	for(int i_reg = 0; i_reg < (int)sample_pe_regions->size(); i_reg++)
	{
		sample_pe_regions->at(i_reg)->score = 1; // This connection is differential.

		char* cur_line = (char*)(sample_pe_regions->at(i_reg)->data);

		// Parse the line.
		// 1       30001   31001   2       114340001       114341001       11
		char dest_chrom[1000];
		int dest_start;
		int dest_end;
		int n_disc_reads = 0;
		double cur_bin_mapp = 0;
		double cur_read_mapp = 0;
		if(sscanf(cur_line, "%*s %*s %*s %lf %lf %s %d %d %d", &cur_bin_mapp, &cur_read_mapp, dest_chrom, &dest_start, &dest_end, &n_disc_reads) != 6)
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		t_annot_region* connecting_region = get_empty_region();
		connecting_region->chrom = t_string::copy_me_str(dest_chrom);
		connecting_region->start = dest_start;
		connecting_region->end = dest_end;
		connecting_region->strand = '+';
		connecting_region->score = n_disc_reads;

		t_per_bin_connection_info* cur_reg_conn_info = new t_per_bin_connection_info();
		cur_reg_conn_info->connecting_bins = new vector<t_connection_info_entry*>();
		t_connection_info_entry* new_conn_info = new t_connection_info_entry();
		new_conn_info->connecting_region = connecting_region;
		cur_reg_conn_info->connecting_bins->push_back(new_conn_info);
		cur_reg_conn_info->avg_mappability = cur_bin_mapp;
		cur_reg_conn_info->avg_read_mappability = cur_read_mapp;

		delete [] cur_line;
		sample_pe_regions->at(i_reg)->data = cur_reg_conn_info;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d sample regions.\n", (int)sample_pe_regions->size());

	// Load the control file.
	vector<t_annot_region*>* control_pe_regions = load_BED_with_line_information(control_fp);
	if(control_pe_regions == NULL)
	{
		fprintf(stderr, "Could not load %s\n", control_fp);
		exit(0);
	}

	for(int i_reg = 0; i_reg < (int)control_pe_regions->size(); i_reg++)
	{
		control_pe_regions->at(i_reg)->score = 1; // This connection is differential.

		char* cur_line = (char*)(control_pe_regions->at(i_reg)->data);

		// Parse the line.
		// 1       30001   31001   2       114340001       114341001       11
		char dest_chrom[1000];
		int dest_start;
		int dest_end;
		int n_disc_reads = 0;
		double cur_bin_mapp = 0;
		double cur_read_mapp = 0;
		if(sscanf(cur_line, "%*s %*s %*s %lf %lf %s %d %d %d", &cur_bin_mapp, &cur_read_mapp, dest_chrom, &dest_start, &dest_end, &n_disc_reads) != 6)
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		t_annot_region* connecting_region = get_empty_region();
		connecting_region->chrom = t_string::copy_me_str(dest_chrom);
		connecting_region->start = dest_start;
		connecting_region->end = dest_end;
		connecting_region->strand = '+';
		connecting_region->score = n_disc_reads;

		t_per_bin_connection_info* cur_reg_conn_info = new t_per_bin_connection_info();
		cur_reg_conn_info->connecting_bins = new vector<t_connection_info_entry*>();
		t_connection_info_entry* new_conn_info = new t_connection_info_entry();
		new_conn_info->connecting_region = connecting_region;
		cur_reg_conn_info->connecting_bins->push_back(new_conn_info);
		cur_reg_conn_info->avg_mappability = cur_bin_mapp;
		cur_reg_conn_info->avg_read_mappability = cur_read_mapp;

		delete [] cur_line;
		control_pe_regions->at(i_reg)->data = cur_reg_conn_info;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d control regions.\n", (int)control_pe_regions->size());

	// Following identifies the regions that do not overlap either the first or the second region in the connection formed by the discordant pair.
	vector<t_annot_region*>* overlaps = intersect_annot_regions(sample_pe_regions, control_pe_regions, true);
	fprintf(stderr, "Found %d overlaps.\n", (int)overlaps->size());
	vector<t_annot_region*>* overlapping_sample_pe_regions = new vector<t_annot_region*>();
	for(int i_int = 0; i_int < (int)overlaps->size(); i_int++)
	{
		t_intersect_info* ol_info = (t_intersect_info*)(overlaps->at(i_int)->data);
		overlapping_sample_pe_regions->push_back(ol_info->src_reg);

		// Check if this connection matches in the second region.
		t_annot_region* sample_second_reg = ((t_per_bin_connection_info*)(ol_info->src_reg->data))->connecting_bins->at(0)->connecting_region;
		t_annot_region* control_second_reg = ((t_per_bin_connection_info*)(ol_info->dest_reg->data))->connecting_bins->at(0)->connecting_region;

		if(t_string::compare_strings(sample_second_reg->chrom, control_second_reg->chrom))
		{
			int ol_start = MAX(control_second_reg->start, sample_second_reg->start);
			int ol_end = MIN(control_second_reg->end, sample_second_reg->end);
			if(ol_end + 2000 > ol_start - 2000)
			{
				// there is overlap, set this region to overlapping.
				ol_info->src_reg->score = 0; // This connection is not differential.
			}
			else
			{
				// Does not overlap in the second region here.
			}
		}
		else
		{
			// Does not overlap in the second region here.
		}
	} // i_int loop.

	//FILE* f_diff_disc_reads = open_f("differential_pe_reads.txt", "w");
	FILE* f_diff_disc_reads = open_f(op_fp, "w");
	for(int i_reg = 0; i_reg < (int)(sample_pe_regions->size()); i_reg++)
	{
		// Is this a differential pe?
		if(sample_pe_regions->at(i_reg)->score == 1)
		{
			t_annot_region* second_region = ((t_per_bin_connection_info*)(sample_pe_regions->at(i_reg)->data))->connecting_bins->at(0)->connecting_region;
			t_annot_region* first_region = sample_pe_regions->at(i_reg);

			fprintf(f_diff_disc_reads, "%s\t%d\t%d\t%lf\t%lf\t%s\t%d\t%d\t%d\n", first_region->chrom, translate_coord(first_region->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(first_region->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				((t_per_bin_connection_info*)(sample_pe_regions->at(i_reg)->data))->avg_mappability,
				((t_per_bin_connection_info*)(sample_pe_regions->at(i_reg)->data))->avg_read_mappability,
				second_region->chrom, second_region->start, second_region->end, second_region->score);
		}
	} // i_diff loop.
	fclose(f_diff_disc_reads);
}

void sort_read_lines_per_id_in_place(vector<char*>* cur_chr_read_lines)
{
	// Replace tab/space with an end-of-string.
	for(int i_read = 0; i_read < (int)(cur_chr_read_lines->size()); i_read++)
	{
		int l_cur_line = t_string::string_length(cur_chr_read_lines->at(i_read));
		bool set_eos = false;
		for(int i = 0; i < l_cur_line; i++)
		{
			if(cur_chr_read_lines->at(i_read)[i] == ' ' || 
				cur_chr_read_lines->at(i_read)[i] == '\t')
			{
				cur_chr_read_lines->at(i_read)[i] = 0;
				set_eos = true;
				break;
			}
		} // i loop.

		if(!set_eos)
		{
			fprintf(stderr, "Did not find a delimited for %s\n", cur_chr_read_lines->at(i_read));
			exit(0);
		}
	} // i_read loop.

	// Sort the read lines.
	fprintf(stderr, "Starting in-place sort\n");
	sort(cur_chr_read_lines->begin(), cur_chr_read_lines->end(), sort_read_lines);
	fprintf(stderr, "Done in-place sort\n");

	// Replace the tabs positions back.
	for(int i_read = 0; i_read < (int)(cur_chr_read_lines->size()); i_read++)
	{
		int i = 0;
		while(1)
		{
			// Repalce the end of string with a tab.
			if(cur_chr_read_lines->at(i_read)[i] == 0)
			{
				cur_chr_read_lines->at(i_read)[i] = '\t';
				break;
			}

			i++;
		} // i loop.
	} // i_read loop.
}

vector<char*>* sort_read_lines_per_id(vector<char*>* cur_chr_read_lines)
{
	vector<t_read_line_w_id*>* read_line_entries = new vector<t_read_line_w_id*>();
	for(int i_read = 0; i_read < (int)(cur_chr_read_lines->size()); i_read++)
	{
		t_read_line_w_id* cur_line_w_id = new t_read_line_w_id();
		char cur_id[1000];
		sscanf(cur_chr_read_lines->at(i_read), "%s", cur_id);
		cur_line_w_id->id = t_string::copy_me_str(cur_id);
		cur_line_w_id->read_line = cur_chr_read_lines->at(i_read);

		read_line_entries->push_back(cur_line_w_id);
	} // i_read loop.

	sort(read_line_entries->begin(), read_line_entries->end(), sort_read_line_entries_per_id);

	vector<char*>* sorted_read_lines = new vector<char*>();
	for(int i_read = 0; i_read < (int)(read_line_entries->size()); i_read++)
	{
		sorted_read_lines->push_back(read_line_entries->at(i_read)->read_line);

		delete [] read_line_entries->at(i_read)->id;
		delete read_line_entries->at(i_read);
	} // i_read loop.
	delete read_line_entries;

	return(sorted_read_lines);
}

vector<char*>* sort_bucket_read_lines(char* bucket_fp)
{
	// Load the reads.
	vector<char*>* bucket_read_lines = buffer_file(bucket_fp);	
	vector<int>* read_starts = new vector<int>();
	vector<t_read_line_sorting_info*>* sorting_info_list = new vector<t_read_line_sorting_info*>();
	for(int i_read = 0; i_read < (int)bucket_read_lines->size(); i_read++)
	{
		int cur_read_start = 0;
		sscanf(bucket_read_lines->at(i_read), "%*s %*s %d", &cur_read_start);

		t_read_line_sorting_info* cur_line_info = new t_read_line_sorting_info();
		cur_line_info->start = cur_read_start;
		cur_line_info->read_line = bucket_read_lines->at(i_read);

		sorting_info_list->push_back(cur_line_info);
	} // i_read loop.

	sort(sorting_info_list->begin(), sorting_info_list->end(), sort_read_line_info);
	vector<char*>* sorted_bucket_read_lines = new vector<char*>();

	for(int i_read = 0; i_read < (int)sorting_info_list->size(); i_read++)
	{
		sorted_bucket_read_lines->push_back(sorting_info_list->at(i_read)->read_line);

		delete sorting_info_list->at(i_read);
	} // i_read loop.

	delete(sorting_info_list);
	delete(read_starts);

	delete(bucket_read_lines);

	return(sorted_bucket_read_lines);
}

// Following is for sorting the mapped reads offline.
bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2)
{
	return(info1->start < info2->start);
}

void load_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads)
{
	// Use the piped input.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = open_f(sam_fp, "r");
	}

	char cur_seq_id[1000];
	char cur_seq_nucs[1000];
	char phred_quality_str[1000];
	while(1)
	{
		char* cur_sam_line = getline(f_sam);
		if(cur_sam_line == NULL)
		{
			break;
		}

		// Parse the sam read line.
		if(sscanf(cur_sam_line, "%s %*d %*s %*d %*d %*s %*s %*d %*d %s %s", cur_seq_id, cur_seq_nucs, phred_quality_str) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_sam_line);
			exit(0);
		}

		// Add the new sequenced read.
		t_sequenced_read* new_sequenced_read = new t_sequenced_read();
		new_sequenced_read->id = t_string::copy_me_str(cur_seq_id);
		new_sequenced_read->nucs = t_string::copy_me_str(cur_seq_nucs);
		new_sequenced_read->mapping_info = NULL;
		new_sequenced_read->quality_str = t_string::copy_me_str(phred_quality_str);
		sequenced_reads->push_back(new_sequenced_read);

		if(sequenced_reads->size() % 100000 == 0)
		{
			fprintf(stderr, "Processing %ld. read.          \r", sequenced_reads->size());
		}
	} // file reading loop.

	// If the sam file is not stdin, close it.
	if(strcmp(sam_fp, "stdin") != 0)
	{
		fclose(f_sam);
	}
}

// FASTQ LOADING INTERFACE:
// Following are for fastq loading from mapped read files.
void load_sequenced_reads_per_fastq(char* fastq_fp, vector<t_sequenced_read*>* sequenced_reads)
{
	// Use the piped input.
	FILE* f_fastq = NULL;
	if(strcmp(fastq_fp, "stdin") == 0)
	{
		f_fastq = stdin;
	}
	else
	{
		f_fastq = open_f(fastq_fp, "r");
	}

	while(1)
	{
		char* cur_seq_id = getline(f_fastq);
		if(cur_seq_id == NULL)
		{
			break;
		}

		char* cur_seq_nucs = getline(f_fastq);
		if(cur_seq_nucs == NULL)
		{
			fprintf(stderr, "Could not read the sequence for %s.\n", cur_seq_id);
			exit(0);
		}

		char* cur_seq_qual_id = getline(f_fastq);
		if(cur_seq_qual_id == NULL)
		{
			fprintf(stderr, "Could not read the quality id for %s.\n", cur_seq_id);
			exit(0);
		}

		char* cur_seq_qual_str = getline(f_fastq);
		if(cur_seq_qual_str == NULL)
		{
			fprintf(stderr, "Could not read the quality str for %s.\n", cur_seq_id);
			exit(0);
		}

		// Add the new sequenced read.
		t_sequenced_read* new_sequenced_read = new t_sequenced_read();
		new_sequenced_read->id = cur_seq_id;
		new_sequenced_read->nucs = cur_seq_nucs;
		new_sequenced_read->mapping_info = NULL;
		new_sequenced_read->quality_str = cur_seq_qual_str;
		sequenced_reads->push_back(new_sequenced_read);

		if(sequenced_reads->size() % 100000 == 0)
		{
			fprintf(stderr, "Processing %ld. read.          \r", sequenced_reads->size());
		}
	} // file reading loop.

	if(strcmp(fastq_fp, "stdin") != 0)
	{
		fclose(f_fastq);
	}
}

void get_per_strand_read_stats(vector<t_annot_region*>* regions,
	char* preprocessed_reads_dir)
{
	t_restr_annot_region_list* restructured_regions = restructure_annot_regions(regions);

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);

	for(int i_reg_chr = 0; i_reg_chr < (int)restructured_regions->chr_ids->size(); i_reg_chr++)
	{
		// Load the fragments for the current chromosome.
		int i_read_chr = t_string::get_i_str(chr_ids, restructured_regions->chr_ids->at(i_reg_chr));

		if(i_read_chr < (int)chr_ids->size())
		{
			char mapped_reads_fp[1000];
			//sprintf(mapped_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_read_chr));
			if (!set_preprocessed_read_file_path_per_dir_chr(preprocessed_reads_dir, chr_ids->at(i_read_chr), mapped_reads_fp))
			{
				fprintf(stderr, "Could not find preprocessed reads file for %s in %s.\n", preprocessed_reads_dir, chr_ids->at(i_read_chr));
				exit(0);
			}

			// Load the fragments for the current chromosome.
			vector<t_mapped_fragment*>* fore_strand_fragments = new vector<t_mapped_fragment*>();
			vector<t_mapped_fragment*>* rev_strand_fragments = new vector<t_mapped_fragment*>();
			load_fragments(mapped_reads_fp, 
				fore_strand_fragments, rev_strand_fragments, 
				0);

			// Sort the fragments before setting the 3p ends.
			sort(fore_strand_fragments->begin(), fore_strand_fragments->end(), sort_mapped_fragments);
			sort(rev_strand_fragments->begin(), rev_strand_fragments->end(), sort_mapped_fragments);

            vector<int>* cur_chr_fore_fragment_3p_posns = new vector<int>();
            for(int i = 0; i < (int)fore_strand_fragments->size(); i++)
            {
				cur_chr_fore_fragment_3p_posns->push_back(fore_strand_fragments->at(i)->base_index + fore_strand_fragments->at(i)->sequenced_fragment_length - 1);
            } // i loop.

            vector<int>* cur_chr_rev_fragment_3p_posns = new vector<int>();
            for(int i = 0; i < (int)rev_strand_fragments->size(); i++)
            {
				//cur_chr_rev_fragment_3p_posns->push_back(rev_strand_fragments->at(i)->base_index);
				cur_chr_rev_fragment_3p_posns->push_back(rev_strand_fragments->at(i)->base_index + rev_strand_fragments->at(i)->sequenced_fragment_length - 1);
				//rev_strand_fragments->at(i)->base_index -= (rev_strand_fragments->at(i)->sequenced_fragment_length - 1);
            } // i loop.

			// Go over all the regions in this chromosome.
			vector<t_annot_region*>* cur_chr_regions = restructured_regions->regions_per_chrom[i_reg_chr];
			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
			{
				double cur_reg_n_mapped_fore_fragments = 0.0;
				double cur_reg_n_mapped_fore_nucs = 0.0;
				get_read_statistics_per_region(fore_strand_fragments, 
													cur_chr_fore_fragment_3p_posns,
													cur_chr_regions->at(i_reg), 
													cur_reg_n_mapped_fore_fragments, 
													cur_reg_n_mapped_fore_nucs);

				double cur_reg_n_mapped_rev_fragments = 0.0;
				double cur_reg_n_mapped_rev_nucs = 0.0;
				get_read_statistics_per_region(rev_strand_fragments, 
													cur_chr_rev_fragment_3p_posns, 
													cur_chr_regions->at(i_reg), 
													cur_reg_n_mapped_rev_fragments, 
													cur_reg_n_mapped_rev_nucs);

				// Set the read stats to the region.
				double* cur_cnts = new double[2];
				cur_cnts[0] = cur_reg_n_mapped_fore_nucs;
				cur_cnts[1] = cur_reg_n_mapped_rev_nucs;
				cur_chr_regions->at(i_reg)->data = cur_cnts;
			} // i_reg loop.

			// Free memory for the current chromosome.
			delete_fragments(fore_strand_fragments);
			delete_fragments(rev_strand_fragments);
			delete(cur_chr_fore_fragment_3p_posns);
			delete(cur_chr_rev_fragment_3p_posns);
		} // read chr check.
	} // i_reg_chr loop.
}

void delete_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads)
{
	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		delete [] sequenced_reads->at(i_r)->id;
		delete [] sequenced_reads->at(i_r)->nucs;
		delete [] sequenced_reads->at(i_r)->quality_str;

		delete sequenced_reads->at(i_r);
	} // i_r loop.

	delete sequenced_reads;
}

void dump_fastq(vector<t_sequenced_read*>* sequenced_reads, char* op_fastq_fp)
{
	FILE* f_op = open_f(op_fastq_fp, "w");

	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		fprintf(f_op, "@%s\n%s\n+%s\n%s\n", 
			&(sequenced_reads->at(i_r)->id[1]), 
			sequenced_reads->at(i_r)->nucs,
			&(sequenced_reads->at(i_r)->id[1]), 
			sequenced_reads->at(i_r)->quality_str);
	} // i_r loop.

	fclose(f_op);
}

void subsample_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads, 
	vector<t_sequenced_read*>* subsampled_sequenced_reads, 
	int n_subsample_size)
{
	if((int)sequenced_reads->size() < n_subsample_size)
	{
		fprintf(stderr, "The size of the sequenced reads is smaller than the requested subsample size: %ld, %d\n", sequenced_reads->size(), n_subsample_size);
		exit(0);
	}
	else
	{
		fprintf(stderr, "Sampling %d reads from the list of %ld reads.\n", n_subsample_size, sequenced_reads->size());
	}

	char* chosen_flag = new char[sequenced_reads->size()+1];
	memset(chosen_flag, 0, sequenced_reads->size());

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	double choosing_threshold = (double)n_subsample_size / (int)sequenced_reads->size();
	fprintf(stderr, "Sampling threshold is %lf\n", choosing_threshold);

	// Make sure the requested subsample size is generated.
	while((int)subsampled_sequenced_reads->size() < n_subsample_size)
	{
		// Go over all the reads.
		for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
		{
			if((int)subsampled_sequenced_reads->size() >= n_subsample_size)
			{
				break;
			}

			// Did we choose this read before?
			if(!chosen_flag[i_r])
			{
				double cur_rand = rng->random_double_ran3();
				if(cur_rand < choosing_threshold)
				{
					chosen_flag[i_r] = true;
					subsampled_sequenced_reads->push_back(sequenced_reads->at(i_r));

					if(subsampled_sequenced_reads->size() % 1000000 == 0)
					{
						fprintf(stderr, "Sampled %ld. read.            \r", subsampled_sequenced_reads->size());
					}
				}
				else
				{
					// This read will not be chosen.
				}
			}
			else
			{
				// Do not try to add the read again.
			}
		} // i_r loop.
	} // subsample size check.

	delete [] chosen_flag;
	delete(rng);
}

void dump_phred_quality_distribution(vector<t_sequenced_read*>* sequenced_reads, char* op_fp)
{
	vector<vector<char>*>* distributions_per_read_posn = new vector<vector<char>*>();
	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		if(i_r % 100000 == 0)
		{
			fprintf(stderr, "Processing %d. read.           \r", i_r);
		}

		// Process the reads.
		int l_cur_read = strlen(sequenced_reads->at(i_r)->nucs);

		// Add the entries if it is necessary to extend the posns.
		if(l_cur_read >= (int)distributions_per_read_posn->size())
		{
			for(int i = (int)distributions_per_read_posn->size(); i < l_cur_read; i++)
			{
				fprintf(stderr, "Adding %d. position to the distribution.\n", i);
				distributions_per_read_posn->push_back(new vector<char>());
			} // i loop.
		}

		// Process all the entries.
		for(int i = 0; i < l_cur_read; i++)
		{
			distributions_per_read_posn->at(i)->push_back(sequenced_reads->at(i_r)->quality_str[i]);
		} // i loop.
	} // i_r loop.

	// Dump the mean and variance for the qualities.
	FILE* f_quals = open_f(op_fp, "w");
	for(int i = 0; i < (int)distributions_per_read_posn->size(); i++)
	{
		vector<char>* cur_qual_sample = distributions_per_read_posn->at(i);

		double qual_sum = 0.0;
		for(int i_q = 0; i_q < (int)cur_qual_sample->size(); i_q++)
		{
			qual_sum += cur_qual_sample->at(i_q);
		} // i_q loop.

		double mean = qual_sum / (int)cur_qual_sample->size();
		double var = 0.0;
		for(int i_q = 0; i_q < (int)cur_qual_sample->size(); i_q++)
		{
			var += (cur_qual_sample->at(i_q) - mean) * (cur_qual_sample->at(i_q) - mean);
		} // i_q loop.

		var /= (cur_qual_sample->size()-1);

		fprintf(f_quals, "%d\t%lf\t%lf\n", i, mean, var);
	} // i loop.
	fclose(f_quals);
}

void load_mapped_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_ELAND(char* eland_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_tagAlign(char* tagalign_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_bowtie(char* bowtie_fp, vector<t_sequenced_read*>* sequenced_reads);
// ABOVE IS THE FASTQ LOADING INTERFACE:

//struct t_PE_read_info
//{
//	char* read_id;
//	char* chr_id;
//	char*
//};

void preprocess_PE_fragments_per_pos_sorted_SAM_file(char* mrf_fp, char* parsed_reads_op_dir, void (preprocess_mapped_read_line)(char* cur_line,
													char* read_id,
													char* chrom,
													int& chr_index, int& sequenced_length,
													char& strand_char,
													char* mapping_quality_str),
													int max_l_fragment,
													bool dump_read_id)
{
	FILE* f_mrf = open_f(mrf_fp, "r");

	if (f_mrf == NULL)
	{
		fprintf(stderr, "Could not open %s\n", mrf_fp);
		return;
	}

	int n_reads = 0;

	int MAX_N_READS_IN_BUFFER = 2000;

	vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
	vector<char*>* frag_fps = new vector<char*>();

	vector<char*>* buffered_read_ids = new vector<char*>();
	vector<char*>* buffered_read_chrs = new vector<char*>();
	vector<int>* buffered_read_posns = new vector<int>();
	vector<int>* buffered_read_flags = new vector<int>();
	vector<char>* buffered_read_strands = new vector<char>();

	// Check chromosome id's list file.
	char chr_ids_fp[100000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);

	vector<char*>* chr_ids = NULL;
	if(check_file(chr_ids_fp))
	{
		chr_ids = buffer_file(chr_ids_fp);

		fprintf(stderr, "Found chromosome id's @ %s, pooling.\n", chr_ids_fp);

		// Open the files for appending.
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			char new_fn[1000];
			sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i_chr));
			
			if(check_file(new_fn))
			{
				fprintf(stderr, "Opening %s for pooling.\n", new_fn);
				frag_f_ptrs->push_back(open_f(new_fn, "a"));
				frag_fps->push_back(t_string::copy_me_str(new_fn));
			}
			else
			{
				sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chr_ids->at(i_chr));
				if (!check_file(new_fn))
				{
					fprintf(stderr, "Could not open %s\n", new_fn);
					open_f(chr_ids_fp, "w");
					exit(0);
				}

				fprintf(stderr, "Opening %s for pooling.\n", new_fn);
				frag_f_ptrs->push_back(open_f(new_fn, "a"));
				frag_fps->push_back(t_string::copy_me_str(new_fn));
			}
		} // i_chr loop.
	}
	else
	{
		// The chromosomes will be added now.
		chr_ids = new vector<char*>();
	}

	if (frag_fps->size() != frag_f_ptrs->size())
	{
		fprintf(stderr, "Sanity check failed: Number of files do not match to pointers.\n");
		exit(0);
	}

	while(1)
	{
		char* cur_line = getline(f_mrf);
		if(cur_line == NULL)
		{
			break;
		}

		n_reads++;

		if (n_reads % 10000 == 0)
		{
			fprintf(stderr, "@%d. read (%d reads in buffer)              \r", n_reads, buffered_read_ids->size());
		}
		
		if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
		{
			fprintf(stderr, "Read: %s (%d reads in buffer)\n", cur_line, buffered_read_ids->size());
			getc(stdin);
		}

		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char; 
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
									read_id,
									chrom, 
									chr_index, sequenced_length, 
									strand_char, 
									mapping_quality_str);

		// Make sure that the line is valid.
		if(chr_index >= 1 &&
			chrom[0] != 0)
		{
			char flag_str[100];
			if (sscanf(cur_line, "%*[^\t] %[^\t]", flag_str) != 1)
			{
				fprintf(stderr, "Could not parse flag string from line: %s\n", cur_line);
				exit(0);
			}

			int flag = atoi(flag_str);

			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if(i_chr == (int)chr_ids->size())
			{
				// Add the chromosome id.
				chr_ids->push_back(t_string::copy_me_str(chrom));
				i_chr = t_string::get_i_str(chr_ids, chrom);

				char new_fn[10000];
				sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chrom);

				// Does the file exist? If so, use the file, do not overwrite.
				frag_f_ptrs->push_back(open_f(new_fn, "w"));
				frag_fps->push_back(t_string::copy_me_str(new_fn));

				fprintf(stderr, "Added %s\n", chrom);
			}

			FILE* cur_frag_file = frag_f_ptrs->at(i_chr);

			if(cur_frag_file == NULL)
			{
				//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
			}
			else
			{
				// Find the current read within the reads in buffer.
				bool found_mate = false;
				for (int i_buff_read = 0; i_buff_read < buffered_read_ids->size(); i_buff_read++)
				{
					if (t_string::compare_strings(buffered_read_chrs->at(i_buff_read), chrom) &&
						buffered_read_posns->at(i_buff_read) > chr_index)
					{
						fprintf(stderr, "The reads seems to be not sorted.\n");
						exit(0);
					}

					// Compare buffered and current read ids and match the distance, and the chromosome id's
					if (t_string::compare_strings(buffered_read_ids->at(i_buff_read), read_id) &&
						t_string::compare_strings(buffered_read_chrs->at(i_buff_read), chrom) &&
						fabs((double)(chr_index + sequenced_length - buffered_read_posns->at(i_buff_read))) < max_l_fragment )
					{
						// Found the matching fragment id.
						char frag_strand_char = strand_char;
						if (buffered_read_flags->at(i_buff_read) & 0x40)
						{
							frag_strand_char = buffered_read_strands->at(i_buff_read);
						}
						
						int start_pos = buffered_read_posns->at(i_buff_read);
						int end_pos = chr_index + sequenced_length;

						// Save: 
						fprintf(cur_frag_file, "%dM %c %d\n", end_pos - start_pos, frag_strand_char, start_pos);

						found_mate = true;

						if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
						{
							fprintf(stderr, "Found matching left mate: %s @ %d\n", buffered_read_ids->at(i_buff_read), buffered_read_posns->at(i_buff_read));
							getc(stdin);
						}

						// Remove the left mapping mate information.
						delete[] buffered_read_ids->at(i_buff_read);
						delete[] buffered_read_chrs->at(i_buff_read);
						buffered_read_ids->erase(buffered_read_ids->begin() + i_buff_read);
						buffered_read_chrs->erase(buffered_read_chrs->begin() + i_buff_read);
						buffered_read_posns->erase(buffered_read_posns->begin() + i_buff_read);
						buffered_read_strands->erase(buffered_read_strands->begin() + i_buff_read);
						buffered_read_flags->erase(buffered_read_flags->begin() + i_buff_read);

						break;
					}
				} // i_buff_read loop.

				if (!found_mate)
				{
					buffered_read_ids->push_back(t_string::copy_me_str(read_id));
					buffered_read_chrs->push_back(t_string::copy_me_str(chrom));
					buffered_read_posns->push_back(chr_index);
					buffered_read_strands->push_back(strand_char);
					buffered_read_flags->push_back(flag);
				}

				// Remove the reads that could not be matched.
				int i_buff_read_2_remove_per_very_long_frag = -1;
				int end_pos = chr_index + sequenced_length;
				for (int i_buff_read = 0; i_buff_read < buffered_read_ids->size(); i_buff_read++)
				{
					// If the chromosomes do not match or the maximum fragment length is not holding, remove.
					if (!t_string::compare_strings(buffered_read_chrs->at(i_buff_read), chrom) ||
						(end_pos - buffered_read_posns->at(i_buff_read) > max_l_fragment))
					{
						delete[] buffered_read_ids->at(i_buff_read);
						delete[] buffered_read_chrs->at(i_buff_read);
						i_buff_read_2_remove_per_very_long_frag = i_buff_read;
					}
					else
					{
						break;
					}
				} // i_buff_read loop.

				// Report the number of reads that are discarded at the end of the chromosome.
				if (buffered_read_chrs->size() > 0)
				{
					if (!t_string::compare_strings(buffered_read_chrs->at(0), chrom))
					{
						fprintf(stderr, "Discarding %d reads after moving to chromosome %s\n", i_buff_read_2_remove_per_very_long_frag+1, chrom);
					}
				}

				if (i_buff_read_2_remove_per_very_long_frag >= 0)
				{
					buffered_read_ids->erase(buffered_read_ids->begin(), buffered_read_ids->begin() + i_buff_read_2_remove_per_very_long_frag + 1);
					buffered_read_chrs->erase(buffered_read_chrs->begin(), buffered_read_chrs->begin() + i_buff_read_2_remove_per_very_long_frag + 1);
					buffered_read_posns->erase(buffered_read_posns->begin(), buffered_read_posns->begin() + i_buff_read_2_remove_per_very_long_frag + 1);
					buffered_read_strands->erase(buffered_read_strands->begin(), buffered_read_strands->begin() + i_buff_read_2_remove_per_very_long_frag + 1);
					buffered_read_flags->erase(buffered_read_flags->begin(), buffered_read_flags->begin() + i_buff_read_2_remove_per_very_long_frag + 1);
				}

				/*if(dump_read_id)
				{
					fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
				}
				else
				{
					fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
				}*/
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete [] cur_line;
	} // file reading loop.

	// (Re)Dump the chromosome id list.
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for(int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);

	// Close fragment file pointers.
	for(int i_f = 0; i_f < (int)frag_f_ptrs->size(); i_f++)
	{
		close_f(frag_f_ptrs->at(i_f), frag_fps->at(i_f));

		delete[] frag_fps->at(i_f);
	} // i_f loop.

	delete frag_fps;
	delete frag_f_ptrs;

	// Unload/close the mapped read file.
	close_f(f_mrf, mrf_fp);
}

// This function open the last file; optimal for sorted data.
void preprocess_mapped_reads_file_single_file_buffering(char* mrf_fp, char* parsed_reads_op_dir, 
	void (preprocess_mapped_read_line)(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str),
	bool dump_read_id)
{
	FILE* f_mrf = open_f(mrf_fp, "r");
	if (f_mrf == NULL)
	{
		fprintf(stderr, "Could not open %s\n", mrf_fp);
		return;
	}

	//char cur_line[100000];
	int n_frags = 0;
	//int n_total_frags = 0;

	FILE* cur_preprocessed_read_f_ptr = NULL;
	char* cur_preprocessed_read_fp = NULL;
	char* cur_preprocessed_chr = NULL;
	vector<char*>* chr_ids = new vector<char*>();

	while (1)
	{
		char* cur_line = getline(f_mrf);
		if (cur_line == NULL)
		{
			break;
		}

		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char;
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);

		// Make sure that the line is valid.
		if (chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if (i_chr == (int)chr_ids->size())
			{
				// Add the chromosome id.
				chr_ids->push_back(t_string::copy_me_str(chrom));

				fprintf(stderr, "Added %s\n", chrom);
			}

			// Select the file pointer.
			if (cur_preprocessed_chr == NULL ||
				(cur_preprocessed_chr != NULL && !t_string::compare_strings(cur_preprocessed_chr, chrom)))
			{
				// Copy the current chromosome.
				if (cur_preprocessed_chr != NULL)
				{
					delete[] cur_preprocessed_chr;
				}

				cur_preprocessed_chr = t_string::copy_me_str(chrom);

				// If the last file is open, close it.
				if (cur_preprocessed_read_f_ptr != NULL)
				{
					// Close the file.
					fprintf(stderr, "Closing %s\n", cur_preprocessed_read_fp);
					close_f(cur_preprocessed_read_f_ptr, cur_preprocessed_read_fp);
				}

				// Re-open the file.
				if (cur_preprocessed_read_fp == NULL)
				{
					cur_preprocessed_read_fp = new char[1000];
				}

				sprintf(cur_preprocessed_read_fp, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chrom);
				cur_preprocessed_read_f_ptr = open_f(cur_preprocessed_read_fp, "a");
			}

			FILE* cur_frag_file = cur_preprocessed_read_f_ptr;

			if (cur_frag_file == NULL)
			{
				//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
			}
			else
			{
				if (dump_read_id)
				{
					fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
				}
				else
				{
					fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
				}
				n_frags++;
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete[] cur_line;
	} // file reading loop.

	// Unload/close the mapped read file.
	close_f(f_mrf, mrf_fp);

	if (cur_preprocessed_read_f_ptr != NULL)
	{
		close_f(cur_preprocessed_read_f_ptr, cur_preprocessed_read_fp);
	}

	// (Re)Dump the chromosome id list.
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for (int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);
}

// Generic preprocessing function for mapped read files.
void preprocess_mapped_reads_file(char* mrf_fp, 
	char* parsed_reads_op_dir, 
	vector<char*>* preset_chr_ids,
	void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	bool dump_read_id)
{
    FILE* f_mrf = open_f(mrf_fp, "r");
	if(f_mrf == NULL)
	{
		fprintf(stderr, "Could not open %s\n", mrf_fp);
		return;
	}

    //char cur_line[100000];
    int n_frags = 0;
    //int n_total_frags = 0;

    vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
	vector<char*>* frag_fps = new vector<char*>();

	// Check chromosome id's list file.
	char chr_ids_fp[100000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);

	vector<char*>* chr_ids = NULL;
	if (preset_chr_ids == NULL)
	{
		if (check_file(chr_ids_fp))
		{
			chr_ids = buffer_file(chr_ids_fp);

			fprintf(stderr, "Found chromosome id's @ %s, pooling.\n", chr_ids_fp);

			// Open the files for appending.
			for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
			{
				char new_fn[1000];
				sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chr_ids->at(i_chr));
				frag_f_ptrs->push_back(open_f(new_fn, "a"));
				frag_fps->push_back(t_string::copy_me_str(new_fn));
			} // i_chr loop.
		}
		else
		{
			// The chromosomes will be added now.
			chr_ids = new vector<char*>();
		}
	} // preset_chr_ids check.
	else
	{
		fprintf(stderr, "Using preset %d chromosomes.\n", preset_chr_ids->size());

		// Preset chromosomes exist.
		chr_ids = new vector<char*>();
		for (int i_chr = 0; i_chr < preset_chr_ids->size(); i_chr++)
		{
			char* cur_chr_id = t_string::copy_me_str(preset_chr_ids->at(i_chr));
			normalize_chr_id(cur_chr_id);

			chr_ids->push_back(cur_chr_id);

			char new_fn[1000];
			sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, cur_chr_id);

			fprintf(stderr, "Opening %s for writing.\n", new_fn);
			frag_f_ptrs->push_back(open_f(new_fn, "a"));
			frag_fps->push_back(t_string::copy_me_str(new_fn));
		} // i_chr loop.
	} // preset chromosome check.

	if (frag_fps->size() != frag_f_ptrs->size())
	{
		fprintf(stderr, "Sanity check failed: Number of files do not match to pointers.\n");
		exit(0);
	}

	while(1)
	{
		char* cur_line = getline(f_mrf);
		if(cur_line == NULL)
		{
			break;
		}
		
		n_frags++;

		if (n_frags % (1000 * 1000) == 0)
		{
			fprintf(stderr, "@ %d reads.                  \r", n_frags);
		}

		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char; 
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
									read_id,
									chrom, 
									chr_index, sequenced_length, 
									strand_char, 
									mapping_quality_str);

		// Make sure that the line is valid.
		if(chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			bool chromosome_preset_confirm = true;
			if(i_chr == (int)chr_ids->size())
			{
				// If there were no preset chromosome id's, update the list.
				if (preset_chr_ids == NULL)
				{
					// Add the chromosome id.
					chr_ids->push_back(t_string::copy_me_str(chrom));
					i_chr = t_string::get_i_str(chr_ids, chrom);

					char new_fn[10000];
					sprintf(new_fn, "%s/%s_mapped_reads.txt.gz", parsed_reads_op_dir, chrom);

					// Re-open this file from scratch; even if it was there, it will be overwritten.
					frag_f_ptrs->push_back(open_f(new_fn, "w"));
					frag_fps->push_back(t_string::copy_me_str(new_fn));

					fprintf(stderr, "Added %s\n", chrom);
				}
				else
				{
					// Skip adding this chromosome; it is not in the preset list.
					chromosome_preset_confirm = false;

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Skipping processing %s               \n", cur_line);
					}
				}
			}

			if (chromosome_preset_confirm)
			{
				FILE* cur_frag_file = frag_f_ptrs->at(i_chr);

				if (cur_frag_file == NULL)
				{
					//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
				}
				else
				{
					if (dump_read_id)
					{
						fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
					}
					else
					{
						fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
					}
				}
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete [] cur_line;
	} // file reading loop.

	// (Re)Dump the chromosome id list.
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for(int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);

	// Close fragment file pointers.
	for(int i_f = 0; i_f < (int)frag_f_ptrs->size(); i_f++)
	{
		close_f(frag_f_ptrs->at(i_f), frag_fps->at(i_f));

		// Compress, if it is already not compressed.
		if (!t_string::ends_with(frag_fps->at(i_f), "gz"))
		{
			char comp_frag_fp[1000];
			sprintf(comp_frag_fp, "%s.gz", frag_fps->at(i_f));
			fprintf(stderr, "Compressing to %s\n", comp_frag_fp);
			compressFile(frag_fps->at(i_f), comp_frag_fp);
			delete_file(frag_fps->at(i_f));
		}

		delete[] frag_fps->at(i_f);
	} // i_f loop.

	delete frag_fps;
	delete frag_f_ptrs;

	// Unload/close the mapped read file.	
	close_f(f_mrf, mrf_fp);
}

void preprocess_mapped_PE_SAM_file(char* pe_sam_fp,
	int min_mapping_qual, 
	char* first_reads_dir, char* last_reads_dir)
{
    // Divide SAM output with respect to chromosomes.
    FILE* f_pe_sam = NULL;
	if(strcmp(pe_sam_fp, "stdin") == 0)
	{
		f_pe_sam = stdin;
	}
	else
	{
		f_pe_sam = open_f(pe_sam_fp, "r");
	}

	if(f_pe_sam == NULL)
	{
		fprintf(stderr, "Could not open %s\n", pe_sam_fp);
		return;
	}

    //char cur_line[100000];
    int n_frags = 0;
    //int n_total_frags = 0;

    vector<FILE*>* frag_1_f_ptrs = new vector<FILE*>();
	vector<FILE*>* frag_2_f_ptrs = new vector<FILE*>();

	// Check chromosome id's list file.
	char chr_ids_fp[100000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", first_reads_dir);

	vector<char*>* chr_ids = NULL;
	//if(check_file(chr_ids_fp))
	//{
	//	chr_ids = buffer_file(chr_ids_fp);

	//	fprintf(stderr, "Found chromosome id's @ %s, pooling.\n", chr_ids_fp);

	//	// Open the files for appending.
	//	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	//	{
	//		char new_fn[1000];
	//		sprintf(new_fn, "%s/%s_mapped_reads_1.txt", parsed_reads_op_dir, chr_ids->at(i_chr));
	//		frag_1_f_ptrs->push_back(open_f(new_fn, "a"));
	//		sprintf(new_fn, "%s/%s_mapped_reads_2.txt", parsed_reads_op_dir, chr_ids->at(i_chr));
	//		frag_2_f_ptrs->push_back(open_f(new_fn, "a"));
	//	} // i_chr loop.
	//}
	//else
	//{
	// The chromosomes will be added now.
	chr_ids = new vector<char*>();
	//}

	// Start reading the file.
	while(1)
	{
		char* cur_line = getline(f_pe_sam);

		if(cur_line == NULL)
		{
			break;
		}
		
		// Load the mapping info based on the file type.
		char cur_read_id[1000];
		char chrom[1000];
		int chr_index = 0;
		int sequenced_length = 0;
		char strand_char = 0; 
		char cigar_str[20000];
		int mapq;
		bool first_segment_in_template;
		bool last_segment_in_template;
		preprocess_PE_SAM_read_line(cur_line,
									cur_read_id, 
									chrom, 
									first_segment_in_template, 
									last_segment_in_template, 
									chr_index, sequenced_length, 
									strand_char,
									mapq,
									cigar_str);

		// Make sure that the line is valid.
		if((min_mapping_qual < 0 || min_mapping_qual < mapq) &&
			chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the 
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if(i_chr == (int)chr_ids->size())
			{
				// Add the chromosome id.
				chr_ids->push_back(t_string::copy_me_str(chrom));
				i_chr = t_string::get_i_str(chr_ids, chrom);

				char new_fn_1[10000];
				sprintf(new_fn_1, "%s/%s_mapped_reads.txt", first_reads_dir, chrom);
				char new_fn_2[10000];
				sprintf(new_fn_2, "%s/%s_mapped_reads.txt", last_reads_dir, chrom);

				// Does the file exist? If so, use the file, do not overwrite.
				if(check_file(new_fn_1))
				{					
					if(!check_file(new_fn_2))
					{
						fprintf(stderr, "Could not find %s.\n", new_fn_2);
						exit(0);
					}

					fprintf(stderr, "Found %s, %s, concatting %s\n", new_fn_1, new_fn_2, chrom);

					fprintf(stderr, "Concatting to %s\n", new_fn_1);
					frag_1_f_ptrs->push_back(open_f(new_fn_1, "a"));
					frag_2_f_ptrs->push_back(open_f(new_fn_2, "a"));
				}
				else
				{
					fprintf(stderr, "Added %s\n", chrom);
					frag_1_f_ptrs->push_back(open_f(new_fn_1, "w"));
					frag_2_f_ptrs->push_back(open_f(new_fn_2, "w"));
				}				
			}

			FILE* cur_frag_file = NULL;
			if(first_segment_in_template)
			{
				cur_frag_file = frag_1_f_ptrs->at(i_chr);
			}
			else if(last_segment_in_template)
			{
				cur_frag_file = frag_2_f_ptrs->at(i_chr);
			}

			if(cur_frag_file == NULL)
			{
					//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
			}
			else
			{
				fprintf(cur_frag_file, "%s %s %c %d\n", cur_read_id, cigar_str, strand_char, chr_index);				
				n_frags++;
			}
		} // check if the line corresponds to a valid mapped nucleotide.
		else
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Skipping : %s\n", cur_line);
			getc(stdin);
}
		}

		delete [] cur_line;
	} // file reading loop.

	// (Re)Dump the chromosome id list.
	if(check_file(chr_ids_fp))
	{
		vector<char*>* existing_chr_ids = buffer_file(chr_ids_fp);

		// Append to the existing chromosome id's.
		FILE* f_chrs = open_f(chr_ids_fp, "a");
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			int i_str = t_string::get_i_str(existing_chr_ids, chr_ids->at(i_chr));
			if(i_str == (int)existing_chr_ids->size())
			{
				fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
			}
		} // i_chr loop.
		fclose(f_chrs);	
	}
	else
	{
		FILE* f_chrs = open_f(chr_ids_fp, "w");
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
		} // i_chr loop.

		fclose(f_chrs);
	}

	// Close fragment file pointers.
	for(int i_f = 0; i_f < (int)frag_1_f_ptrs->size(); i_f++)
	{
		fclose(frag_1_f_ptrs->at(i_f));
	}

	for(int i_f = 0; i_f < (int)frag_2_f_ptrs->size(); i_f++)
	{
		fclose(frag_2_f_ptrs->at(i_f));
	}

	// Unload/close the mapped read file.
	if(strcmp(pe_sam_fp, "stdin") == 0)
	{
		fclose(f_pe_sam);
	}
}

void count_mapped_reads_per_file(char* mrf_fp, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	double n_total_reads)
{
    // Divide SAM output with respect to chromosomes.
	FILE* f_mrf = open_f(mrf_fp, "r");

	if(f_mrf == NULL)
	{
		fprintf(stderr, "mapped read file pointer and file buffer are both NULL for %s\n", mrf_fp);
		return;
	}

    //char cur_line[2000];
    //int n_frags = 0;
    //int n_total_frags = 0;
	n_total_reads = 0;

	fprintf(stderr, "Counting the mapped reads from %s\n", mrf_fp);
	while(1)
	{
		//char* cur_line = getline(f_mrf);
		char* cur_line = getline(f_mrf);

		if(cur_line == NULL)
		{
			break;
		}
		
		// Load the mapping info based on the file type.
		//char chrom[1000];
		//int chr_index;
		//int sequenced_length;
		//char strand_char; 
		//char mapping_quality_str[1000];
		//preprocess_mapped_read_line(cur_line, 
		//							chrom, 
		//							chr_index, sequenced_length, 
		//							strand_char, 
		//							mapping_quality_str);
		char phred_quality_str[1000];
		int flag;
		if(cur_line[0] == '@')
		{
		}
		else
		{
			flag = 0;
			//if(sscanf(cur_line, "%*s %*d %*s %*d %*d %*s %*s %*d %*d %*s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 1)
			if(sscanf(cur_line, "%*s %d %*s %*d %*d %*s %*s %*d %*d %*s %s", &flag, phred_quality_str) == 2)
			{
				// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
				//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
				//chr_index = translate_coord(chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

				// Check the flag and determine the strand.
				/*strand_char = 'F';
				if(flag & 0x10)
				{
					strand_char = 'R';
				}*/

				// Sanity check. Is this fragment mapped?
				if(flag & 0x04)
				{
					// The read is not mapping.
					//chrom[0] = 0;
				}
				else
				{
					n_total_reads++;
				}
			}
			else
			{
				// Could not parse the line.
			}
		}

		if( (((int)n_total_reads) % 1000000) == 0)
		{
			fprintf(stderr, "At %lf. read.\n", n_total_reads);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "%lf reads.\n", n_total_reads);

	// Unload/close the mapped read file.
	close_f(f_mrf, mrf_fp);
}

void preprocess_tagAlign_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	int chr_start_index;
	int chr_end_index;
	char strand_sign;

	if(sscanf(cur_line, "%s %d %d %*s %*d %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
	{
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		//chr_start_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		//chr_end_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		chr_start_index = translate_coord(chr_start_index, TAGALIGN_COORDS::start_base, CODEBASE_COORDS::start_base);
		chr_end_index = translate_coord(chr_end_index, TAGALIGN_COORDS::end_base, CODEBASE_COORDS::end_base);

		// Set quality to all matches.
		sprintf(cigar_str, "%dM", chr_end_index-chr_start_index+1);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = chr_end_index-chr_start_index+1;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_preprocessed_LH_GFF3_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int _chr_index;
	char strand;

	//if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 6)
	// X       14705460        35M     +
	if(sscanf(cur_line, "%s %d %s %c", chrom, &_chr_index, mapping_quality_str, &strand) == 4)
	{
		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		chr_index = _chr_index;

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand == '-')
		{
			strand_char = 'R';
		}

			chr_index = _chr_index;
			sequenced_length = 0;
	}
	else
	{
		chrom[0] = 0;
	}
}

// TODO: Make an option to load only a certain phed partition from the file, so that we conserve memory.
unsigned short*** load_partitioned_compressed_pileups(char* cur_comp_allele_fp, int& n_partitions, int& l_pileup)
{
	FILE* f_comp = open_f(cur_comp_allele_fp, "rb");

	// Load the partitions.
	int n_read_partitions = 0;
	fread(&n_read_partitions, sizeof(int), 1, f_comp);
	n_partitions = n_read_partitions;

	// Load the partitions.
	int l_read_sig = 0;
	fread(&l_read_sig, sizeof(int), 1, f_comp);
	int l_sig = l_read_sig;

	l_pileup = l_sig;

	fprintf(stderr, "Loading pileup of length %d over %d phred partitions\n", l_sig, n_partitions);

	unsigned short*** loaded_pileup = new unsigned short**[n_partitions];

	for(int part_i = 0; part_i < n_partitions; part_i++)
	{
		loaded_pileup[part_i] = allocate_pileup(l_sig);
	} // part_i loop.

	// Go over all the positions.
	int i = 1;
	while(i <= l_sig)
	{
		if(i % 1000000 == 0)
		{
			fprintf(stderr, "Loading %d            \r", i);
		}

		// Read the existence flag.
		unsigned char existence_flag = 0;
		fread(&existence_flag, sizeof(unsigned char), 1, f_comp);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
		fprintf(stderr, "Loading %d. position, existence flag: 0x%x\n", i, existence_flag);

		// Check an RLE case.
		if(existence_flag == 0xFF)
		{
			unsigned int l_RLE = 0;
			fread(&l_RLE, sizeof(unsigned int), 1, f_comp);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Loading RLE of length %d @ %d\n", l_RLE, i);

			// When we add this, we move upto the location where a 0-run ends.
			i += l_RLE;
		} // RLE check.
		else
		{
			// Read the first byte, parse the data range and which alleles exist.
			unsigned char drange_flag = ((existence_flag & (1 << 5)) >> 5);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Loading val @ %d: drange_flag: %d\n", i, drange_flag);

			for(int part_i = 0; part_i < n_partitions; part_i++)
			{
				for(int allele_i = 0; allele_i < 5; allele_i++)
				{
					// Do we have an entry here?
					if((existence_flag & (1 << allele_i)) > 0)
					{
						if(drange_flag == 1)
						{	
							unsigned short current_count_short = 0;
							fread(&(current_count_short), sizeof(unsigned short), 1, f_comp);
							loaded_pileup[part_i][allele_i][i] = current_count_short;
						}
						else
						{
							unsigned char current_count_char = 0;
							fread(&(current_count_char), sizeof(unsigned char), 1, f_comp);
							loaded_pileup[part_i][allele_i][i] = (unsigned short)current_count_char;
						}

						//fprintf(stderr, "%d, ", loaded_pileup[part_i][allele_i][i]);
					} // count check for current position.
				} // allele_i loop.
			} // part_i loop.

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "\n");

			// Update the position.
			i++;
		} // non-RLE check.
	} // i loop.
	close_f(f_comp, cur_comp_allele_fp);

	// Free memory.
	return(loaded_pileup);
}

unsigned short** allocate_pileup(int l_sig)
{
	unsigned short** loaded_pileup = new unsigned short*[5];
	for (int i_allele = 0; i_allele < 5; i_allele++)
	{
		loaded_pileup[i_allele] = new unsigned short[l_sig + 2];
		memset(loaded_pileup[i_allele], 0, sizeof(unsigned short) * (l_sig + 2));
	} // i_allele loop.

	return loaded_pileup;
}

void delete_pileup(unsigned short** loaded_pileup)
{
	//unsigned short** loaded_pileup = new unsigned short*[5];
	for (int i_allele = 0; i_allele < 5; i_allele++)
	{
		delete [] loaded_pileup[i_allele];
	} // i_allele loop.

	delete[] loaded_pileup;
}

// We are assuming a consolidated indel callset at this point with filtered and homopolymer/repeat resolved.
void compute_single_cell_allelic_stats_per_10X_SAM(char* per_cell_barcodes_fp, char* SAM_fp, char* chr_info_list_fp, char* genome_seq_dir, char* candidate_snvs_bed_fp, char* candidate_indels_bed_fp, char* op_fp)
{
	fprintf(stderr, "Counting the ref/alt alleles per cell from %s using variants in %s (snv) and %s (indel).\n", SAM_fp, candidate_snvs_bed_fp, candidate_indels_bed_fp);
	fprintf(stderr, "Insertion coordinates are: [chrom] [Nucleotide coordinate to the left of insert] [Nucleotide coordinates to the right of insert] [Ref/alt string]\n\
Deletion coordinates are: [chrom] [Leftmost deleted nucleotide coordinate] [Rightmost deleted nucleotide coordinate] [Ref/alt string]\n");

	vector<char*>* per_cell_barcodes = buffer_file(per_cell_barcodes_fp);
	sort(per_cell_barcodes->begin(), per_cell_barcodes->end(), t_string::sort_strings_per_prefix);
	fprintf(stderr, "Loaded and sorted %d cell barcodes.\n", per_cell_barcodes->size());

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chr_info_list_fp, chr_ids, chr_lengths);

	// All variant regions should be like : 
	vector<t_annot_region*>* snv_regs = new vector<t_annot_region*>();
	if (check_file(candidate_snvs_bed_fp))
	{
		snv_regs = load_BED(candidate_snvs_bed_fp);
	} // snv file existence check.
	fprintf(stderr, "Loaded %d SNVs.\n", snv_regs->size());

	// Load the indels.
	vector<t_annot_region*>* indel_regs = new vector<t_annot_region*>();
	if (check_file(candidate_indels_bed_fp))
	{
		indel_regs = load_BED(candidate_indels_bed_fp);
	} // indel file existence check.
	fprintf(stderr, "Loaded %d indels.\n", indel_regs->size());

	// Pool variants.
	vector<t_annot_region*>* pooled_var_regs = new vector<t_annot_region*>();
	pooled_var_regs->insert(pooled_var_regs->end(), snv_regs->begin(), snv_regs->end());
	pooled_var_regs->insert(pooled_var_regs->end(), indel_regs->begin(), indel_regs->end());
	fprintf(stderr, "Setting up counts for %d pooled variants.\n", pooled_var_regs->size());

	// Set the information for each variant.
	vector<t_annot_region*>* ins_regs = new vector<t_annot_region*>();
	vector<t_annot_region*>* del_regs = new vector<t_annot_region*>();
	for (int i_var = 0;
		i_var < pooled_var_regs->size();
		i_var++)
	{
		t_string_tokens* refalt_toks = t_string::tokenize_by_chars(pooled_var_regs->at(i_var)->name, " ");
		if (refalt_toks->size() != 2)
		{
			fprintf(stderr, "The refalt allele string in variant region name is not as expected @ %s(%d): %s\n", __FILE__, __LINE__, pooled_var_regs->at(i_var)->name);
			exit(0);
		}

		const char* ref_all_str = refalt_toks->at(0)->str();
		const char* alt_all_str = refalt_toks->at(1)->str();

		t_variant_per_cell_stats* var_stats = new t_variant_per_cell_stats();

		if (ref_all_str[0] == '.')
		{
			var_stats->var_type = VAR_TYPE_INSERTION;
			ins_regs->push_back(pooled_var_regs->at(i_var));
		}
		else if (alt_all_str[0] == '.')
		{
			var_stats->var_type = VAR_TYPE_DELETION;
			del_regs->push_back(pooled_var_regs->at(i_var));
		}
		else if(t_string::string_length(ref_all_str) == 1 &&
				t_string::string_length(alt_all_str) == 1)
		{
			var_stats->var_type = VAR_TYPE_SNV;

			if (pooled_var_regs->at(i_var)->start != pooled_var_regs->at(i_var)->end)
			{
				fprintf(stderr, "SNV start is not equal to end: %s:%d-%d: %s\n", 
						pooled_var_regs->at(i_var)->chrom,
						pooled_var_regs->at(i_var)->start,
						pooled_var_regs->at(i_var)->end,
						pooled_var_regs->at(i_var)->name);
			}
		}
		else
		{
			fprintf(stderr, "Could not determine variant type from refalt string: %s: %s, %s\n", pooled_var_regs->at(i_var)->name, ref_all_str, alt_all_str);
			exit(0);
		}

		var_stats->ref_allele = t_string::copy_me_str(ref_all_str);
		var_stats->alt_allele = t_string::copy_me_str(alt_all_str);
		var_stats->alt_allele_cnt_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(var_stats->alt_allele_cnt_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));
		var_stats->ref_allele_cnt_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(var_stats->ref_allele_cnt_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		// Parse the ref/alt alleles.
		pooled_var_regs->at(i_var)->data = var_stats;

		t_string::clean_tokens(refalt_toks);
	} // i_snv loop.

	fprintf(stderr, "Setup %d SNV, %d INS, %d DEL regions.\n", snv_regs->size(), ins_regs->size(), del_regs->size());

	// We do not need to load chromosomes any more since the alt/ref allele information is stored in the variant files.
	//fprintf(stderr, "Loading chromosome sequences.\n");
	//char** per_chrom_seq = new char*[chr_ids->size() + 2];
	//int* per_chrom_lengths = new int[chr_ids->size() + 2];
	//for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	//{
	//	int l_loaded_seq = 0;
	//	char cur_chrom_seq_fp[1000];
	//	sprintf(cur_chrom_seq_fp, "%s/%s.bin", genome_seq_dir, chr_ids->at(i_chr));
	//	if (!check_file(cur_chrom_seq_fp))
	//	{
	//		sprintf(cur_chrom_seq_fp, "%s/%s.bin.gz", genome_seq_dir, chr_ids->at(i_chr));

	//		if (!check_file(cur_chrom_seq_fp))
	//		{
	//			fprintf(stderr, "Could not find the sequence file %s\n", cur_chrom_seq_fp);
	//			exit(0);
	//		}
	//	}

	//	per_chrom_seq[i_chr] = load_binary_sequence_file(cur_chrom_seq_fp, l_loaded_seq);
	//	per_chrom_lengths[i_chr] = l_loaded_seq;
	//	fprintf(stderr, "Loaded %s (%d, %d)                   \r", chr_ids->at(i_chr), l_loaded_seq, chr_lengths->at(i_chr));
	//} // i_chr loop.

	fprintf(stderr, "Restructuring the pooled variants into %d chromosomes.\n", chr_ids->size());
	t_restr_annot_region_list* restr_pooled_regs = restructure_annot_regions(pooled_var_regs, chr_ids);

	// Set the sorting information that is necessary for overlapping reads.
	for (int i_chr = 0; i_chr < restr_pooled_regs->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(restr_pooled_regs->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	fprintf(stderr, "Done, Reading %s ..\n", SAM_fp);

	// Enter file reading loop.
	char cur_fragment[1000];
	char phred_quality_str[100000];

	unsigned long long int n_total_processed_reads = 0;
	unsigned long long int n_snv_containing_reads = 0;
	unsigned long long int n_ins_containing_reads = 0;
	unsigned long long int n_del_containing_reads = 0;
	unsigned long long int n_unmapped_reads = 0;

	unsigned long long int n_matching_matching_nucs = 0;
	unsigned long long int n_total_matching_nucs = 0;

	unsigned long long int n_reads_w_non_matched_bcs = 0;
	unsigned long long int n_no_barcode_reads = 0;

	FILE* f_sam = open_f(SAM_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_sam);
		if (cur_line == NULL)
		{
			break;
		}

		// This is the global parameter, update this before checking for thread outof/which check.
		n_total_processed_reads++;

		// Report only with the lead thread.
		if (n_total_processed_reads % (1000 * 1000) == 0)
		{
			fprintf(stderr, "Processing %llu. read: SNV: %llu; INS: %llu; DEL: %llu; Unmapped: %llu; NOBC: %llu; UNMATCHEDBC: %llu; (%llu/%llu matching nucs).               \r",
					n_total_processed_reads,
					n_snv_containing_reads,
					n_ins_containing_reads,
					n_del_containing_reads,
					n_unmapped_reads,
					n_no_barcode_reads,
					n_reads_w_non_matched_bcs,
					n_matching_matching_nucs,
					n_total_matching_nucs);
		}

		// If this is a comment line, skip it.
		if (cur_line[0] == '@')
		{
			delete[] cur_line;
			continue;
		}

		char read_id[1000];
		char chrom[100];
		int chr_index;
		int sequenced_length;
		char strand_char;
		char cigar_str[100];

		preprocess_SAM_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			cigar_str);

		// Copy the cigar string.
		char* mapping_map_str = cigar_str;

		// Make sure that the line is valid.
		if (chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int chr_i = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if (chr_i == (int)chr_ids->size())
			{
				// This read is on a chromosome we do not care about. Count it then throw away?
			}
			else
			{
				// Parse and search the cellular barcode and locate the cell, first.
				int cell_i = -1;
				char barcode_buffer[1000];
				bool found_CB_flag = get_10X_cell_barcode_per_SAM_read(cur_line, barcode_buffer);
				if (!found_CB_flag)
				{
					n_no_barcode_reads++;

					// Could not parse the barcode; free memory and continue.
					delete[] cur_line;
					continue;
				}
				else
				{
					int search_cell_i = t_string::fast_search_string_per_prefix(barcode_buffer, per_cell_barcodes, 0, per_cell_barcodes->size() - 1);
					while (search_cell_i > 0 &&
						(t_string::sort_strings_per_prefix(barcode_buffer, per_cell_barcodes->at(search_cell_i)) ||
							t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
					{
						search_cell_i--;
					} // search_cell_i loop.

					while (search_cell_i < per_cell_barcodes->size() &&
						(t_string::sort_strings_per_prefix(per_cell_barcodes->at(search_cell_i), barcode_buffer) ||
							t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
					{
						if (t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer))
						{
							cell_i = search_cell_i;
							break;
						}
						else
						{
							search_cell_i++;
						}
					} // search_cell_i loop.
				} // barcode check.

				if (cell_i == -1)
				{
					n_reads_w_non_matched_bcs++;
					delete[] cur_line;
					continue;
				}

				// Load the sequence and quality data.
				if (sscanf(cur_line, "%*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]",
					cur_fragment, phred_quality_str) == 2)
				{

				}
				else
				{
					fprintf(stderr, "Could not parse %s\n", cur_line);
					exit(0);
				}

				// Parse the cigar string.
				int i_mapp_map = 0;
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				while (mapping_map_str_valid &&
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
						i_mapp_map,
						is_matching,
						l_cur_entry,
						entry_type_char);

					// Find the variant that overlaps with this block.
					// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
					vector<t_annot_region*>* cur_chr_pooled_var_regs = restr_pooled_regs->regions_per_chrom[chr_i];
					int i_var = locate_posn_region_per_region_starts(chr_index, cur_chr_pooled_var_regs, 0, cur_chr_pooled_var_regs->size() - 1);

					// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
					while (i_var > 0 && cur_chr_pooled_var_regs->at(i_var)->sort_info->cumulative_sorted_end >= chr_index)
					{
						i_var--;
					} // i_reg loop.

					// Now go over all the regions and identify the matching.
					while (i_var < cur_chr_pooled_var_regs->size() &&
							cur_chr_pooled_var_regs->at(i_var)->start <= (chr_index + l_cur_entry - 1))
					{
						int ol_start = MAX(cur_chr_pooled_var_regs->at(i_var)->start, chr_index);
						int ol_end = MIN(cur_chr_pooled_var_regs->at(i_var)->end, chr_index + l_cur_entry - 1);

						// Check the overlap; does this block overlap with the coordinates of this variant?
						if (ol_start <= ol_end)
						{
							t_variant_per_cell_stats* cur_var_info = (t_variant_per_cell_stats*)(cur_chr_pooled_var_regs->at(i_var)->data);

							if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
							{
								// Write the mapping between the block and the ref/alt alleles.
								fprintf(stderr, "Block %c:%d (%d)\nVariant %s:%d-%d; REF:%s, ALT:%s ; ",
									entry_type_char, chr_index, l_cur_entry,
									cur_chr_pooled_var_regs->at(i_var)->chrom, cur_chr_pooled_var_regs->at(i_var)->start, cur_chr_pooled_var_regs->at(i_var)->end, 
									cur_var_info->ref_allele, cur_var_info->alt_allele);

								// Build the indicator string.
								t_string* indicator_str = new t_string();
								indicator_str->sprintf("%s (%s)\n", cur_fragment, mapping_map_str);
								for (int i = 0; i < read_nuc_index; i++)
								{
									indicator_str->concat_char(' ');
								} // i loop.

								  // Move on the current block.
								for (int i = chr_index; i < cur_chr_pooled_var_regs->at(i_var)->start; i++)
								{
									indicator_str->concat_char(' ');
								} // i loop.

								// Write the indicator's position.
								for (int i = cur_chr_pooled_var_regs->at(i_var)->start; i <= cur_chr_pooled_var_regs->at(i_var)->end; i++)
								{
									indicator_str->concat_char('^');
								} // i loop.

								fprintf(stderr, "\n%s\n", indicator_str->str());
								delete(indicator_str);
							} // block-2-refalt allele dump,

							// We process only the matching blocks for the SNVs.
							if (cur_var_info->var_type == VAR_TYPE_SNV)
							{
								if (is_matching)
								{
									n_snv_containing_reads++;

									// Update the ref/alt allele on this block.
									int snv_genome_i = cur_chr_pooled_var_regs->at(i_var)->start;
									int snv_read_nuc_i = read_nuc_index + (snv_genome_i - chr_index);

									int ref_cnt = (nuc_2_num(cur_fragment[snv_read_nuc_i]) == nuc_2_num(cur_var_info->ref_allele[0])) ? (1) : (0);
									int alt_cnt = (nuc_2_num(cur_fragment[snv_read_nuc_i]) == nuc_2_num(cur_var_info->alt_allele[0])) ? (1) : (0);

									// Update this cell's allele contributed by this read.
									if (ref_cnt > 0)
									{
										cur_var_info->ref_allele_cnt_per_cell[cell_i]++;

										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "SNV REF.\n");
											getc(stdin);
										}
									} // ref check.
									else if (alt_cnt > 0)
									{
										cur_var_info->alt_allele_cnt_per_cell[cell_i]++;

										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "SNV ALT.\n");
											getc(stdin);
										}
									} // alt check
									else if (ref_cnt == 0 && alt_cnt == 0)
									{
										// This can happen when there are sequencing errors.
										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "No alt or ref allele for %s:%d (%s, %s) on read:\n%s (%s): \n",
												cur_chr_pooled_var_regs->at(i_var)->chrom, cur_chr_pooled_var_regs->at(i_var)->start,
												cur_var_info->ref_allele, cur_var_info->alt_allele,
												cur_fragment, mapping_map_str);
										}
									} // no alt or ref check.
								} // match block check
								else if(entry_type_char == 'D' ||
										entry_type_char == 'I')
								{
									// If this is a del or ins block, update the ref count.
									if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
									{
										fprintf(stderr, "SNV REF.\n");
										getc(stdin);
									}

									cur_var_info->ref_allele_cnt_per_cell[cell_i]++;
								} // match block check
							} // SNV check.
							else if (cur_var_info->var_type == VAR_TYPE_INSERTION)
							{
								n_ins_containing_reads++;

								char* cur_ins_alt_allele = cur_var_info->alt_allele;
								int l_alt_allele = t_string::string_length(cur_ins_alt_allele);

								// Only insertion blocks can update the alternate allele counts: Following matches the coordinates to the extracted blocks.
								if (entry_type_char == 'I' &&
									l_cur_entry == l_alt_allele &&
									chr_index == cur_chr_pooled_var_regs->at(i_var)->end) // This matches ending coordinate exactly to where we want them; i.e., the insertion is engulfed.
								{
									// Inserts are defined by their insertion position on the genome. We will not do anything with the genome, though.
									int l_frag = t_string::string_length(cur_fragment);

									int ins_allele_match = 1;
									for (int cur_read_i = read_nuc_index;
										cur_read_i <= read_nuc_index + l_cur_entry - 1;
										cur_read_i++)
									{
										// Do a check on the alternate allele for the insertion.
										if ((cur_read_i - read_nuc_index) < l_alt_allele &&
											(cur_read_i < l_frag) &&
											nuc_2_num(cur_ins_alt_allele[cur_read_i - read_nuc_index]) != nuc_2_num(cur_fragment[cur_read_i]))
										{
											ins_allele_match = 0;
											break;
										}
									} // cur_genome_i loop.

									//// Match the read to the alternate allele.
									//int ins_allele_match = 1;
									//int ins_allele_trace_nuc_i = 0;
									//while (ins_allele_trace_nuc_i < l_alt_allele &&
									//		(read_nuc_index + ins_allele_trace_nuc_i) < l_frag)
									//{
									//	if (nuc_2_num(cur_ins_alt_allele[ins_allele_trace_nuc_i]) != nuc_2_num(cur_fragment[read_nuc_index + ins_allele_trace_nuc_i]))
									//	{
									//		ins_allele_match = 0;
									//		break;
									//	}

									//	ins_allele_trace_nuc_i++;
									//} // ins_allele_trace_nuc_i loop.

									// Update the insertion's alt and ref counts based on the match with the alternate allele on the read.
									if (ins_allele_match == 1)
									{
										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "INSERT ALT.\n");
											getc(stdin);
										}

										cur_var_info->alt_allele_cnt_per_cell[cell_i]++;
									}
									if (ins_allele_match == 0)
									{
										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "INSERT REF.\n");
											getc(stdin);
										}

										cur_var_info->ref_allele_cnt_per_cell[cell_i]++;
									}
								} // insertion type block check.
								else if(is_matching)
								{
									if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
									{
										fprintf(stderr, "INSERT REF.\n");
										getc(stdin);
									}

									// This block does not engulfs the variant or it is not an insertion block.
									cur_var_info->ref_allele_cnt_per_cell[cell_i]++;
								}
							} // INS check.
							else if (cur_var_info->var_type == VAR_TYPE_DELETION)
							{
								n_del_containing_reads++;

								// If this is an exact matching deletion at this location, it is alternate allele.
								if (entry_type_char == 'D')
								{
									// Exact match.
									if (chr_index == cur_chr_pooled_var_regs->at(i_var)->start &&
										(chr_index+l_cur_entry-1) == cur_chr_pooled_var_regs->at(i_var)->end)
									{
										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "DELETION ALT (EXACT).\n");
											getc(stdin);
										}

										cur_var_info->alt_allele_cnt_per_cell[cell_i]++;
									}
									else if (ol_start == cur_chr_pooled_var_regs->at(i_var)->start &&
											ol_end == cur_chr_pooled_var_regs->at(i_var)->end)
									{
										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "DELETION ALT (ENGULF).\n");
											getc(stdin);
										}

										// This is not an exact matching deletion.
										cur_var_info->alt_allele_cnt_per_cell[cell_i]++;
									}
									else
									{
										if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
										{
											fprintf(stderr, "DELETION REF.\n");
											getc(stdin);
										}

										// This is not an exact matching deletion.
										cur_var_info->ref_allele_cnt_per_cell[cell_i]++;
									}
								} // del type block check.
								else if(is_matching)
								{
									if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
									{
										fprintf(stderr, "DELETION REF.\n");
										getc(stdin);
									}

									cur_var_info->ref_allele_cnt_per_cell[cell_i]++;
								}
							} // DEL check.
							else
							{
								fprintf(stderr, "We are not supposed to be here @ %s(%d)\n", __FILE__, __LINE__);
								exit(0);
							}
						} // overlap check.

						i_var++;
					} // i_var loop.

					// Update the base for the current entry.
					if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if (check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.
			} // chromosome check.
		} // chr_index check.
		else
		{
			n_unmapped_reads++;
		}

		delete[] cur_line;
	} // sam file reading loop.

	// Close the SAM file.
	close_f(f_sam, SAM_fp);

	// Save the per cell snv/ins/del counts.
	FILE* f_op = open_f(op_fp, "w");

	// Write the header.
	fprintf(f_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
	{
		fprintf(f_op, "\t%s", per_cell_barcodes->at(i_cell));
	} // i_cell loop.
	fprintf(f_op, "\n");
	for (int i_var = 0; i_var < pooled_var_regs->size(); i_var++)
	{
		t_variant_per_cell_stats* cur_var_info = (t_variant_per_cell_stats*)(pooled_var_regs->at(i_var)->data);

		fprintf(f_op, "%s\t%d\t%d\t%s", 
				pooled_var_regs->at(i_var)->chrom,
				translate_coord(pooled_var_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(pooled_var_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				pooled_var_regs->at(i_var)->name);

		for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
		{
			fprintf(f_op, "\t%d %d", 
				cur_var_info->ref_allele_cnt_per_cell[i_cell],
				cur_var_info->alt_allele_cnt_per_cell[i_cell]);
		} // i_cell loop.

		fprintf(f_op, "\n");
	} // i_var loop.
	close_f(f_op, op_fp);
} // compute_single_cell_allelic_stats_per_10X_SAM option.

// This code should be eventually combined with the pileup generation code.
void extract_summarize_indel_containing_read_blocks_per_SAM(char* SAM_fp, char* chrom_info_fp, char* genome_seq_dir, char* op_dir)
{
	fprintf(stderr, "Summarizing and extracting indel supporting reads from %s and saving to %s.\n", SAM_fp, op_dir);

	// Load the chromosome information.
	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chrom_info_fp, chr_ids, chr_lengths);

	vector<t_annot_region*>** per_chr_indel_containing_read_blocks = new vector<t_annot_region*>*[chr_ids->size() + 2];
	char** per_chrom_seq = new char*[chr_ids->size() + 2];
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		per_chr_indel_containing_read_blocks[i_chr] = new vector<t_annot_region*>();

		int l_loaded_seq = 0;
		char cur_chrom_seq_fp[1000];
		sprintf(cur_chrom_seq_fp, "%s/%s.bin", genome_seq_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chrom_seq_fp))
		{
			sprintf(cur_chrom_seq_fp, "%s/%s.bin.gz", genome_seq_dir, chr_ids->at(i_chr));

			if (!check_file(cur_chrom_seq_fp))
			{
				fprintf(stderr, "Could not find the sequence file %s\n", cur_chrom_seq_fp);
				exit(0);
			}
		}

		per_chrom_seq[i_chr] = load_binary_sequence_file(cur_chrom_seq_fp, l_loaded_seq);
		fprintf(stderr, "Loaded %s (%d, %d)                   \r", chr_ids->at(i_chr), chr_lengths->at(i_chr), l_loaded_seq);
	} // i_chr loop.

	// Enter file reading loop.
	char cur_fragment[1000];
	char phred_quality_str[100000];

	unsigned long long n_total_processed_reads = 0;
	unsigned long long n_indel_containing_reads = 0;
	unsigned long long n_unmapped_reads = 0;

	double n_matching_matching_nucs = 0;
	double n_total_matching_nucs = 0;

	FILE* f_sam = open_f(SAM_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_sam);
		if (cur_line == NULL)
		{
			break;
		}

		// This is the global parameter, update this before checking for thread outof/which check.
		n_total_processed_reads++;

		// Report only with the lead thread.
		if (n_total_processed_reads % (1000 * 1000) == 0)
		{
			fprintf(stderr, "@ %llu. read: %llu indel containing; %llu unmapped; Ref match/unmatch million nucs: %.1f/%.1f.               \r",
					n_total_processed_reads,
					n_indel_containing_reads,
					n_unmapped_reads,
					n_matching_matching_nucs / (1000 * 1000),
					n_total_matching_nucs / (1000 * 1000));
		}

		// If this is a comment line, skip it.
		if (cur_line[0] == '@')
		{
			delete[] cur_line;
			continue;
		}

		char read_id[1000];
		char chrom[100];
		int chr_index;
		int sequenced_length;
		char strand_char;
		char cigar_str[100];

		preprocess_SAM_read_line(cur_line,
								read_id,
								chrom,
								chr_index, sequenced_length,
								strand_char,
								cigar_str);

		// Copy the cigar string.
		char* mapping_map_str = cigar_str;

		// Make sure that the line is valid.
		if (chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int chr_i = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if (chr_i == (int)chr_ids->size())
			{
				// This read is on a chromosome we do not care about. Count it then throw away?
			}
			else
			{
				// Load the sequence and quality data.
				if (sscanf(cur_line, "%*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", 
							cur_fragment, phred_quality_str) == 2)
				{

				}
				else
				{
					fprintf(stderr, "Could not parse %s\n", cur_line);
					exit(0);
				}

				char strand_plus_min = '+';
				if (strand_char == 'R')
				{
					strand_plus_min = '-';
				}

				// Parse the cigar string.
				int i_mapp_map = 0;
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				while (mapping_map_str_valid &&
						mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map,
														is_matching,
														l_cur_entry,
														entry_type_char);

					// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
					if (is_matching)
					{
						int cur_read_i = read_nuc_index;
						for (int cur_genome_i = chr_index;
							cur_genome_i <= chr_index + l_cur_entry - 1;
							cur_genome_i++)
						{
							if (cur_genome_i <= chr_lengths->at(chr_i))
							{
								if (nuc_2_num(cur_fragment[cur_read_i]) == nuc_2_num(per_chrom_seq[chr_i][cur_genome_i]))
								{
									n_matching_matching_nucs++;
								}
								else
								{
									n_total_matching_nucs++;
								}
							}
							else
							{
								fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
									chr_ids->at(chr_i), cur_genome_i, chr_lengths->at(chr_i));

								exit(0);
							}

							cur_read_i++;
						} // cur_genome_i loop.
					} // check if this block is matching.
					else if (entry_type_char == 'D')
					{
						n_indel_containing_reads++;

						t_annot_region* new_del_reg = get_empty_region();
						new_del_reg->chrom = t_string::copy_me_str(chrom);
						new_del_reg->start = translate_coord(chr_index, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_del_reg->end = translate_coord(chr_index + l_cur_entry - 1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_del_reg->strand = strand_plus_min;
						new_del_reg->data = NULL;
						new_del_reg->score = VAR_TYPE_DELETION;

						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						int l_ref_allele = (chr_index + l_cur_entry - 1) - chr_index + 1;
						int l_alt_allele = 1;
						char* cur_del_ref_allele = new char[l_ref_allele + 5];
						memset(cur_del_ref_allele, 0, sizeof(char) * (l_ref_allele + 2));
						for (int cur_genome_i = chr_index; cur_genome_i <= chr_index + l_cur_entry - 1; cur_genome_i++)
						{
							if (cur_genome_i <= chr_lengths->at(chr_i))
							{
								cur_del_ref_allele[cur_genome_i - chr_index] = per_chrom_seq[chr_i][cur_genome_i];
							}
							else
							{
								fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
										chr_ids->at(chr_i), cur_genome_i, chr_lengths->at(chr_i));

								exit(0);
							}
						} // cur_genome_i loop.

						// Set the ref allele to the data for deletions.
						new_del_reg->data = cur_del_ref_allele;

						if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
						{
							fprintf(stderr, "DELETION @ %d (Ref. %s): %s\n", chr_index, cur_del_ref_allele, cur_line);
							getc(stdin);
						}

						// Add the del block region.
						per_chr_indel_containing_read_blocks[chr_i]->push_back(new_del_reg);
					} // deletion check.
					else if (entry_type_char == 'I')
					{
						n_indel_containing_reads++;

						// This points to the position right before the insertion.
						t_annot_region* new_ins_reg = get_empty_region();
						new_ins_reg->chrom = t_string::copy_me_str(chrom);
						new_ins_reg->start = translate_coord(chr_index-1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_ins_reg->end = translate_coord(chr_index, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_ins_reg->strand = strand_plus_min;
						new_ins_reg->data = NULL;
						new_ins_reg->score = VAR_TYPE_INSERTION;

						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						int l_ref_allele = 1;
						int l_alt_allele = (chr_index + l_cur_entry - 1) - chr_index + 1;
						char* cur_ins_alt_allele = new char[l_alt_allele + 5];
						memset(cur_ins_alt_allele, 0, sizeof(char) * (l_alt_allele + 2));

						// Insertion to the reference: This is included to one position.
						for (int cur_read_i = read_nuc_index;
							cur_read_i <= read_nuc_index + l_cur_entry - 1;
							cur_read_i++)
						{
							cur_ins_alt_allele[cur_read_i - read_nuc_index] = cur_fragment[cur_read_i];
						} // cur_genome_i loop.
						// Set the ref allele to the data for deletions.
						new_ins_reg->data = cur_ins_alt_allele;

						if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
						{
							fprintf(stderr, "INSERTION @ %d (Alt %s): %s\n", chr_index, cur_ins_alt_allele, cur_line);
							getc(stdin);
						}

						// Add the current insertion region to the list of regions.
						per_chr_indel_containing_read_blocks[chr_i]->push_back(new_ins_reg);
					} // insert block check.

					// Update the base for the current entry.
					if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if (check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.
			} // chromosome check.
		} // chr_index check.
		else
		{
			n_unmapped_reads++;
		}

		delete[] cur_line;
	} // sam file reading loop.

	close_f(f_sam, SAM_fp);

	// Sort each chromosome, then save.
	fprintf(stderr, "Finished processing reads, saving to %s\n", op_dir);
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "%s: %d indel blocks.\n", chr_ids->at(i_chr), per_chr_indel_containing_read_blocks[i_chr]->size());
		char cur_indel_blocks_fp[1000];
		sprintf(cur_indel_blocks_fp, "%s/%s_indel_blocks.bed.gz", op_dir, chr_ids->at(i_chr));

		FILE* f_indel_blocks = NULL;
		if (check_file(cur_indel_blocks_fp))
		{
			fprintf(stderr, "Opening %s for pooling.\n", cur_indel_blocks_fp);
			f_indel_blocks = open_f(cur_indel_blocks_fp, "a");
		}
		else
		{
			fprintf(stderr, "Opening %s for saving.\n", cur_indel_blocks_fp);
			f_indel_blocks = open_f(cur_indel_blocks_fp, "w");
		}

		for (int i_block = 0; i_block < per_chr_indel_containing_read_blocks[i_chr]->size(); i_block++)
		{
			char* ref_allele = ".";
			char* alt_allele = ".";

			// Check for insert/delete.
			if (per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->score == VAR_TYPE_INSERTION)
			{
				alt_allele = (char*)(per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->data);
			} // insertion check.
			else if(per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->score == VAR_TYPE_DELETION)
			{
				ref_allele = (char*)(per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->data);
			} // deletion check.
			else
			{
				fprintf(stderr, "Sanity check failed, the block size does not make sense: %s:%d-%d\n", per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->chrom,
						per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->start,
						per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->end);
			}

			fprintf(f_indel_blocks, "%s\t%d\t%d\t%s %s\t%d\t%c\n", 
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->chrom, 
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->start,
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->end,
				ref_allele,
				alt_allele,
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->score,
				per_chr_indel_containing_read_blocks[i_chr]->at(i_block)->strand);
		} // i_block loop.
		close_f(f_indel_blocks, cur_indel_blocks_fp);
	} // i_chr loop.
}

bool sort_regions_coords_first_names_second(t_annot_region* reg1, t_annot_region* reg2)
{
	if (reg1->start != reg2->start)
	{
		return(reg1->start < reg2->start);
	}
	else if (reg1->end != reg2->end)
	{
		return(reg1->end < reg2->end);
	}
	else
	{
		// Both starts and ends are the same; check the alleles.
		return t_string::sort_strings(reg1->name, reg2->name);
	}
}

// Scans and saves insertion/deletion candidates from the pileup file. The issue is insertions are not correctly stored in terms of allele information.
void scan_indels_per_summarized_indel_blocks(char* chr_ids_lengths_fp, char* indel_blocks_dir, char* genome_seq_dir,
												char* postv_pileup_dir, char* negtv_pileup_dir,
												int min_covg_per_indel, int min_alternate_covg_per_indel, double min_alternate_freq_per_indel, 
												int l_indel_scanning_window,
												char* op_fp)
{
	fprintf(stderr, "Scanning Indels with %d bp long vicinity window with a minimum %d total and %d alternate read coverage and %.3f min alternate frequency.\n", 
		l_indel_scanning_window, min_covg_per_indel, min_alternate_covg_per_indel, min_alternate_freq_per_indel);

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_fp, chr_ids, chr_lengths);

	vector<char*>* ref_alt_alleles_per_cur_win = new vector<char*>();
	vector<t_annot_region*>* indel_center_win_block_regs = new vector<t_annot_region*>();
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Scanning indels on %s\n", chr_ids->at(i_chr));

		// Load the pileup.
		int l_postv_covg = 0;
		int l_negtv_covg = 0;

		char cur_chr_pileup0_fp[1000];
		sprintf(cur_chr_pileup0_fp, "%s/%s_allele_counts.bin.gz", postv_pileup_dir, chr_ids->at(i_chr));

		char cur_chr_pileup1_fp[1000];
		sprintf(cur_chr_pileup1_fp, "%s/%s_allele_counts.bin.gz", negtv_pileup_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chr_pileup0_fp) ||
			!check_file(cur_chr_pileup1_fp))
		{
			fprintf(stderr, "Could not find %s, %s skipping.\n", cur_chr_pileup0_fp, cur_chr_pileup1_fp);
			continue;
		}

		int* postv_covg = load_coverage_per_compressed_pileup_file(cur_chr_pileup0_fp, l_postv_covg);
		int* negtv_covg = load_coverage_per_compressed_pileup_file(cur_chr_pileup1_fp, l_negtv_covg);
		fprintf(stderr, "Loaded %d, %d coverages for postv and negtv signals.\n", l_postv_covg, l_negtv_covg);

		// Load the genome sequence, first.
		int l_loaded_seq = 0;
		char cur_chrom_seq_fp[1000];
		sprintf(cur_chrom_seq_fp, "%s/%s.bin", genome_seq_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chrom_seq_fp))
		{
			sprintf(cur_chrom_seq_fp, "%s/%s.bin.gz", genome_seq_dir, chr_ids->at(i_chr));

			if (!check_file(cur_chrom_seq_fp))
			{
				fprintf(stderr, "Could not find the sequence file %s\n", cur_chrom_seq_fp);
				exit(0);
			}
		}

		char* cur_chrom_seq = load_binary_sequence_file(cur_chrom_seq_fp, l_loaded_seq);
		fprintf(stderr, "Loaded %d (%d) nucleotides for chromosome sequence of %s.\n", chr_lengths->at(i_chr), l_loaded_seq, chr_ids->at(i_chr));

		char cur_indel_blocks_fp[1000];
		sprintf(cur_indel_blocks_fp, "%s/%s_indel_blocks.bed.gz", indel_blocks_dir, chr_ids->at(i_chr));

		// Load and sort the regions: The ref/alt alleles are in the name of the regions.
		vector<t_annot_region*>* cur_chrom_indel_blocks = load_BED(cur_indel_blocks_fp);
		sort(cur_chrom_indel_blocks->begin(), cur_chrom_indel_blocks->end(), sort_regions_coords_first_names_second);
		fprintf(stderr, "Loaded and sorted %d indel blocks on %s\n", cur_chrom_indel_blocks->size(), chr_ids->at(i_chr));

		// Slide a window and identify the indels that are matching in ref/alt allelles and have good coverage.
		t_annot_region* last_block = NULL;
		for(int win_center_block_i = 0; win_center_block_i < cur_chrom_indel_blocks->size(); win_center_block_i++)
		{
			if (win_center_block_i % 1000 == 0)
			{
				fprintf(stderr, "@ %d. block          \r", win_center_block_i);
			}

			// Make sure that we did not process a block with the same position and alleles centered at this location.
			if (last_block == NULL)
			{
				last_block = cur_chrom_indel_blocks->at(win_center_block_i);
			}
			else if (last_block->start == cur_chrom_indel_blocks->at(win_center_block_i)->start &&
					t_string::compare_strings(last_block->name, cur_chrom_indel_blocks->at(win_center_block_i)->name))
			{
				continue;
			}

			last_block = cur_chrom_indel_blocks->at(win_center_block_i);

			// Reset the buffer of blocks.
			int n_indel_supporting_blocks_in_cur_window = 0;
			int n_pos_supporting_blocks_in_cur_window = 0;
			int n_neg_supporting_blocks_in_cur_window = 0;
			ref_alt_alleles_per_cur_win->clear();

			int start_block_i = win_center_block_i;
			while (start_block_i > 0)
			{
				if (t_string::compare_strings(cur_chrom_indel_blocks->at(start_block_i)->name, cur_chrom_indel_blocks->at(win_center_block_i)->name))
				{
					// The block region's name contains the ref/alt allale strings.
					ref_alt_alleles_per_cur_win->push_back(cur_chrom_indel_blocks->at(start_block_i)->name);

					if (cur_chrom_indel_blocks->at(start_block_i)->strand == '+')
					{
						n_pos_supporting_blocks_in_cur_window++;
					}
					else
					{
						n_neg_supporting_blocks_in_cur_window++;
					}
				}

				if ((cur_chrom_indel_blocks->at(win_center_block_i)->start - cur_chrom_indel_blocks->at(start_block_i-1)->start) < l_indel_scanning_window)
				{
					start_block_i--;
				}
				else
				{
					break;
				}
			} // start_block_i loop.

			int end_block_i = win_center_block_i;
			while (end_block_i+1 < cur_chrom_indel_blocks->size())
			{
				// The block region's name contains the ref/alt allele strings.
				if (end_block_i != win_center_block_i &&
					t_string::compare_strings(cur_chrom_indel_blocks->at(end_block_i)->name, cur_chrom_indel_blocks->at(win_center_block_i)->name))
				{
					ref_alt_alleles_per_cur_win->push_back(cur_chrom_indel_blocks->at(end_block_i)->name);

					if (cur_chrom_indel_blocks->at(end_block_i)->strand == '+')
					{
						n_pos_supporting_blocks_in_cur_window++;
					}
					else
					{
						n_neg_supporting_blocks_in_cur_window++;
					}
				}

				if ((cur_chrom_indel_blocks->at(end_block_i+1)->start - cur_chrom_indel_blocks->at(win_center_block_i)->start) < l_indel_scanning_window)
				{
					end_block_i++;
				}
				else
				{
					break;
				}
			} // end_block_i loop.

			if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
			{
				fprintf(stderr, "Center block: %s:%d (%s)::[%d-%d] @ [%d-%d]\n", 
						cur_chrom_indel_blocks->at(win_center_block_i)->chrom, cur_chrom_indel_blocks->at(win_center_block_i)->start, cur_chrom_indel_blocks->at(win_center_block_i)->name, 
						cur_chrom_indel_blocks->at(start_block_i)->start, cur_chrom_indel_blocks->at(end_block_i)->start,
						start_block_i, end_block_i);
				getc(stdin);
			}

			n_indel_supporting_blocks_in_cur_window = ref_alt_alleles_per_cur_win->size();

			// If there is enough blocks, i.e., reads, sort the alleles then check covg.
			if (n_indel_supporting_blocks_in_cur_window > min_alternate_covg_per_indel)
			{
				// Get the covg.
				int mid_posn = (cur_chrom_indel_blocks->at(win_center_block_i)->start + cur_chrom_indel_blocks->at(win_center_block_i)->end) / 2;
				int total_covg = postv_covg[mid_posn] + negtv_covg[mid_posn];

				// Check the alternate AF.
				double approximate_alternate_covg = (double)n_indel_supporting_blocks_in_cur_window / (double)total_covg;
				if (approximate_alternate_covg > min_alternate_freq_per_indel)
				{
					// Only add this if it matches to the 
					t_annot_region* indel_supporting_block_reg = duplicate_region(cur_chrom_indel_blocks->at(win_center_block_i));
					indel_supporting_block_reg->name = t_string::copy_me_str(cur_chrom_indel_blocks->at(win_center_block_i)->name);
					indel_supporting_block_reg->score = n_indel_supporting_blocks_in_cur_window;

					// Write the covg statistics for this indel.
					int* per_strand_counts = new int[10];
					per_strand_counts[0] = n_pos_supporting_blocks_in_cur_window;
					per_strand_counts[1] = n_neg_supporting_blocks_in_cur_window;
					per_strand_counts[2] = total_covg;
					indel_supporting_block_reg->data = per_strand_counts;

					if (__DUMP_MAPPED_READ_TOOLS_MSGS__)
					{
						fprintf(stderr, "Adding block: %s:%d (%s, %d)\n", indel_supporting_block_reg->chrom, indel_supporting_block_reg->start,
							ref_alt_alleles_per_cur_win->at(0), n_indel_supporting_blocks_in_cur_window);
					}

					indel_center_win_block_regs->push_back(indel_supporting_block_reg);
				} // frequency check.
			} // alternate allele covg check.
		} // win_center_block_i loop.

		delete[] postv_covg;
		delete[] negtv_covg;
		delete[] cur_chrom_seq;
		delete_annot_regions(cur_chrom_indel_blocks);
	} // i_chr loop.

	fprintf(stderr, "Identified %d indels, saving to %s.\n", indel_center_win_block_regs->size(), op_fp);

	FILE* f_op = open_f(op_fp, "w");
	for (int i_indel = 0; i_indel < indel_center_win_block_regs->size(); i_indel++)
	{
		int* per_strand_counts = (int*)(indel_center_win_block_regs->at(i_indel)->data);

		fprintf(f_op, "%s\t%d\t%d\t%s\t%d\t+\t%d\t%d\t%d\n", 
			indel_center_win_block_regs->at(i_indel)->chrom,
			translate_coord(indel_center_win_block_regs->at(i_indel)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(indel_center_win_block_regs->at(i_indel)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			indel_center_win_block_regs->at(i_indel)->name,
			indel_center_win_block_regs->at(i_indel)->score,
			per_strand_counts[0], per_strand_counts[1], per_strand_counts[2]);
	} // i_indel loop.
	close_f(f_op, op_fp);
}

void dump_pileup_SNV_candidates_per_stranded_pileup(char* chr_info_fp, 
													char* pileup_strand_0_dir, char* pileup_strand_1_dir, 
													char* bin_seq_dir, double min_covg, double min_alternate_covg, double min_alternate_freq, double max_strand_imbalance, 
													char* op_fp)
{
	fprintf(stderr, "Calling SNVs from stranded pileups with minimum %.1f total and %.1f alternate read coverage and %.3f min alternate frequency and %.3f max strand imbalance.\n", 
			min_covg, min_alternate_covg, min_alternate_freq, max_strand_imbalance);

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();

	int* nuc_2_num_array = get_per_char_as_nuc_2_num_coding_array();
	load_chromosome_lengths_per_tabbed_file(chr_info_fp, chr_ids, chr_lengths);

	FILE* f_op = open_f(op_fp, "w");
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Procesing %s\n", chr_ids->at(i_chr));

		int l_cur_chr_seq = 0;
		char cur_chr_bin_fp[1000];
		sprintf(cur_chr_bin_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));
		char* cur_chr_seq = load_binary_sequence_file(cur_chr_bin_fp, l_cur_chr_seq);

		char cur_chr_pileup0_fp[1000];
		sprintf(cur_chr_pileup0_fp, "%s/%s_allele_counts.bin.gz", pileup_strand_0_dir, chr_ids->at(i_chr));

		char cur_chr_pileup1_fp[1000];
		sprintf(cur_chr_pileup1_fp, "%s/%s_allele_counts.bin.gz", pileup_strand_1_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chr_pileup0_fp) ||
			!check_file(cur_chr_pileup1_fp))
		{
			fprintf(stderr, "Could not find %s, %s skipping.\n", cur_chr_pileup0_fp, cur_chr_pileup1_fp);
		}
		else
		{
			int l_pileup_per_strand[2];
			unsigned short** per_strand_pileup[2];
			per_strand_pileup[0] = load_compressed_pileups(cur_chr_pileup0_fp, l_pileup_per_strand[0]);
			per_strand_pileup[1] = load_compressed_pileups(cur_chr_pileup1_fp, l_pileup_per_strand[1]);
			fprintf(stderr, "Loaded %d and %d long pileups.\n", l_pileup_per_strand[0], l_pileup_per_strand[1]);

			char per_num_nucs[] = "ACGTN";

			for (int i = 1;
				i <= MIN(l_pileup_per_strand[0], l_pileup_per_strand[1]);
				i++)
			{
				int cur_total_covg = 0;
				int max_alternate_cnt = 0;
				int max_alt_allele_i = -1;
				for (int i_all = 0; i_all < 4; i_all++)
				{
					cur_total_covg += (per_strand_pileup[0][i_all][i] + per_strand_pileup[1][i_all][i]);

					int cur_alt_all_covg = (per_strand_pileup[0][i_all][i] + per_strand_pileup[1][i_all][i]);

					bool is_alternate_allele = (nuc_2_num_array[cur_chr_seq[i]] != i_all);
					if (is_alternate_allele &&
						cur_alt_all_covg > min_alternate_covg &&
						cur_alt_all_covg > max_alternate_cnt)
					{
						max_alternate_cnt = cur_alt_all_covg;
						max_alt_allele_i = i_all;
					}
				} // i_all loop.

				// 
				if (cur_total_covg < min_covg ||
					max_alt_allele_i == -1)
				{
					continue;
				}

				double cur_alternate_freq = (double)(max_alternate_cnt) / (double)(cur_total_covg);
				double cur_strand_imbalance = max_strand_imbalance + 1;

				//if (per_strand_pileup[0][max_alt_allele_i][i] > 0 &&
				//	per_strand_pileup[1][max_alt_allele_i][i] > 0)
				{
					double dbl_strand0_cnt = (double)(per_strand_pileup[0][max_alt_allele_i][i]) + 1;
					double dbl_strand1_cnt = (double)(per_strand_pileup[1][max_alt_allele_i][i]) + 1;

					cur_strand_imbalance = MAX( dbl_strand0_cnt / dbl_strand1_cnt,
												dbl_strand1_cnt / dbl_strand0_cnt);
				}

				if (cur_alternate_freq > min_alternate_freq &&
					cur_strand_imbalance < max_strand_imbalance)
				{
					fprintf(f_op, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\n", chr_ids->at(i_chr), i,
							cur_chr_seq[i], per_num_nucs[max_alt_allele_i],
							max_alternate_cnt, cur_total_covg,
							per_strand_pileup[0][max_alt_allele_i][i],
							per_strand_pileup[1][max_alt_allele_i][i]);
				}
			} // i loop.

			// Free per strand pileup memory.
			delete_pileup(per_strand_pileup[0]);
			delete_pileup(per_strand_pileup[1]);
		} // Check if the pileup is loaded.

		// Free sequence memory.
		delete[] cur_chr_seq;
	} // i_chr loop.

	fclose(f_op);
}


void dump_pileup_SNV_candidates(char* chr_info_fp, char* pileup_dir, char* bin_seq_dir, int min_covg, int min_alternate_covg, double min_alternate_freq, char* op_fp)
{
	fprintf(stderr, "Calling SNVs with minimum %d total and %d alternate read coverage and %.3f min alternate frequency.\n", min_covg, min_alternate_covg, min_alternate_freq);

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();

	int* nuc_2_num_array = get_per_char_as_nuc_2_num_coding_array();
	load_chromosome_lengths_per_tabbed_file(chr_info_fp, chr_ids, chr_lengths);

	FILE* f_op = open_f(op_fp, "w");
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Procesing %s\n", chr_ids->at(i_chr));

		int l_cur_chr_seq = 0;
		char cur_chr_bin_fp[1000];
		sprintf(cur_chr_bin_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));
		char* cur_chr_seq = load_binary_sequence_file(cur_chr_bin_fp, l_cur_chr_seq);

		char cur_chr_pileup_fp[1000];
		sprintf(cur_chr_pileup_fp, "%s/%s_allele_counts.bin.gz", pileup_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chr_pileup_fp))
		{
			fprintf(stderr, "Could not find %s, skipping.\n", cur_chr_pileup_fp);
		}
		else
		{
			int l_pileup = 0;
			unsigned short** cur_chr_pileup = load_compressed_pileups(cur_chr_pileup_fp, l_pileup);
			fprintf(stderr, "Loaded %d long pileup.\n", l_pileup);

			char per_num_nucs[] = "ACGTN";

			for (int i = 1; i <= l_pileup; i++)
			{
				int cur_total_covg = 0;
				int max_alternate_cnt = 0;
				int max_alt_allele_i = -1;
				for (int i_all = 0; i_all < 4; i_all++)
				{
					cur_total_covg += cur_chr_pileup[i_all][i];

					bool is_alternate_allele = (nuc_2_num_array[cur_chr_seq[i]] != i_all);
					if (is_alternate_allele &&
						cur_chr_pileup[i_all][i] > min_alternate_covg &&
						max_alternate_cnt < cur_chr_pileup[i_all][i])
					{
						max_alternate_cnt = cur_chr_pileup[i_all][i];
						max_alt_allele_i = i_all;
					}
				} // i_all loop.

				if (cur_total_covg < min_covg  ||
					max_alt_allele_i == -1)
				{
					continue;
				}

				double cur_alternate_freq = (double)(max_alternate_cnt) / (double)(cur_total_covg);
				if (cur_alternate_freq < min_alternate_freq)
				{
					continue;
				}

				fprintf(f_op, "%s\t%d\t%c\t%c\t%d\t%d\n", chr_ids->at(i_chr), i,
					cur_chr_seq[i], per_num_nucs[max_alt_allele_i],
					max_alternate_cnt, cur_total_covg);
			} // i loop.

			  // Free memory.
			delete_pileup(cur_chr_pileup);
		} // Check if the pileup is loaded.

		// Free sequence memory.
		delete[] cur_chr_seq;
	} // i_chr loop.

	fclose(f_op);
}

unsigned short** load_compressed_pileups(char* cur_comp_allele_fp, int& l_pileup)
{
	FILE* f_comp = open_f(cur_comp_allele_fp, "rb");

	if(f_comp == NULL)
	{
		fprintf(stderr, "Could not open %s\n", cur_comp_allele_fp);
		exit(0);
	}

	int l_read_sig = 0;
	fread(&l_read_sig, sizeof(int), 1, f_comp);
	int l_sig = l_read_sig;

	l_pileup = l_sig;

	fprintf(stderr, "Loading pileup of length %d\n", l_sig);

	unsigned short** loaded_pileup = allocate_pileup(l_sig);

	// Go over all the positions.
	int i = 1;
	while(i <= l_sig)
	{
		// Read the existence flag.
		unsigned char existence_flag = 0;
		int n_read = fread(&existence_flag, sizeof(unsigned char), 1, f_comp);

		// We have reached the EOF before reading whole file.
		if (n_read == 0)
		{
			break;
		}

		// Check an RLE case.
		if(existence_flag == 0xFF)
		{
			unsigned int l_RLE = 0;
			fread(&l_RLE, sizeof(unsigned int), 1, f_comp);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Loading RLE of length %d @ %d\n", l_RLE, i);

			// When we add this, we move upto the location where a 0-run ends.
			i += l_RLE;
		} // RLE check.
		else
		{
			// Read the first byte, parse the data range and which alleles exist.
			unsigned char drange_flag = ((existence_flag & (1 << 5)) >> 5);

			for(int allele_i = 0; allele_i < 5; allele_i++)
			{
				// Do we have an entry here?
				if((existence_flag & (1 << allele_i)) > 0)
				{
					if(drange_flag == 1)
					{	
						unsigned short current_count_short = 0;
						fread(&(current_count_short), sizeof(unsigned short), 1, f_comp);
						loaded_pileup[allele_i][i] = current_count_short;
					}
					else
					{
						unsigned char current_count_char = 0;
						fread(&(current_count_char), sizeof(unsigned char), 1, f_comp);
						loaded_pileup[allele_i][i] = (unsigned short)current_count_char;
					}
				} // count check for current position.
			} // allele_i loop.

			i++;
		}
	} // i loop.

	close_f(f_comp, cur_comp_allele_fp);

	// Free memory.
	return(loaded_pileup);
}

int* load_coverage_per_pileups(unsigned short** pileups, int l_sig)
{
	int* covg_signal = new int[l_sig+2];
	memset(covg_signal, 0, sizeof(int)*(l_sig+2));

	for(int i = 1; i <= l_sig; i++)
	{
		int cur_posn_covg = 0;
		for(int allele_i = 0; allele_i < 5; allele_i++)
		{
			cur_posn_covg += pileups[allele_i][i];
		} // allele_i loop.

		covg_signal[i] = cur_posn_covg;
	} // i loop.

	return(covg_signal);
}

int* load_coverage_per_compressed_partitioned_pileup_file(char* comp_pileup_fp, int& l_sig)
{
	fprintf(stderr, "Loading coverage from compressed pileup file %s.\n", comp_pileup_fp);

	// Load the compressed pileup.
	int n_parts = 0;
	unsigned short*** pileups = load_partitioned_compressed_pileups(comp_pileup_fp, n_parts, l_sig);

	int* covg_profile = new int[l_sig + 5];
	memset(covg_profile, 0, sizeof(int) * (l_sig+2));
	for(int i = 1; i <= l_sig; i++)
	{
		for(int part_i = 0; part_i < n_parts; part_i++)
		{
			for(int all_i = 0; all_i < 5; all_i++)
			{
				covg_profile[i] += pileups[part_i][all_i][i];
			} // all_i loop.
		} // part_i loop.
	} // i loop.

	// Free the pileups, since we do not need them anymore.
	for(int part_i = 0; part_i < n_parts; part_i++)
	{
		for(int all_i = 0; all_i < 5; all_i++)
		{
			delete [] pileups[part_i][all_i];
		} // all_i loop.

		delete [] pileups[part_i];
	} // part_i loop.

	delete [] pileups;

	return(covg_profile);
}

int* load_coverage_per_compressed_pileup_file(char* comp_pileup_fp, int& l_sig)
{
	fprintf(stderr, "Loading coverage from compressed pileup file %s.\n", comp_pileup_fp);

	// Load the compressed pileup.
	unsigned short** pileups = load_compressed_pileups(comp_pileup_fp, l_sig);

	// Get the coverage profile
	int* covg_profile = load_coverage_per_pileups(pileups, l_sig);

	// Free the pileups, since we do not need them anymore.
	for(int all_i = 0; all_i < 5; all_i++)
	{
		delete [] pileups[all_i];
	} // all_i loop.

	delete [] pileups;

	return(covg_profile);
}

void compress_partitioned_nucleotide_pileup_track(unsigned short*** pileup, int n_partitions, int l_sig, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	
	// Write the # of phred partitions
	fwrite(&n_partitions, sizeof(int), 1, f_op);

	// Write the length of pileup
	fwrite(&l_sig, sizeof(int), 1, f_op);	
	int n_4bit_posns = 0;
	int n_8bit_posns = 0;
	int n_12bit_posns = 0;
	int n_14bit_posns = 0;
	int n_16bit_posns = 0;
	int n_0_signal = 0;

	int cur_0_run_start = -1;
	for(int i = 1; i <= l_sig; i++)
	{
		// Dump the current 5 levels: 		
		// Generate the existence flag.
		//if(i % 1000 == 0)
		//{
		//	fprintf(stderr, "Position %d:\n", i);
		//}

		unsigned char existence_flag = 0;
		unsigned char drange_flag = 0;
		unsigned short max_val = 0;
		unsigned int total_val = 0;
		for(int part_i = 0; part_i < n_partitions; part_i++)
		{
			for(int allele_i = 0; allele_i < 5; allele_i++)
			{
				if(pileup[part_i][allele_i][i] > 0)
				{	
					total_val += (unsigned int)(pileup[part_i][allele_i][i]);

					if(max_val < pileup[part_i][allele_i][i])
					{
						max_val = pileup[part_i][allele_i][i];
					}

					existence_flag = existence_flag | (1 << allele_i);
					//fprintf(stderr, "Nuc %d: %d\n", allele_i, pileup[allele_i][i]);
					//getc(stdin);
				} // signal check.
			} // allele_i loop.
		} // part_i loop.

		// Following is an RLE check.
		if(total_val == 0)
		{
			if(cur_0_run_start == -1)
			{
				// Initiate the 0 run.
				cur_0_run_start = i;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
				fprintf(stderr, "Initiated a 0-run @ %d\n", i);
			}
			else
			{
				// We are already in a zero run.
			}
		} // total_val 0 check.
		else
		{
			// Check if we will do a RLE dump.
			// Make sure there was a run of zeros before this position.
			if(cur_0_run_start != -1)
			{
				bool do_RLE_dump = ((i - cur_0_run_start) > 5);

				if(do_RLE_dump)
				{
					//fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

					// Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
					// Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
					unsigned char RLE_indicator = 0xFF;
					fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
					fprintf(stderr, "Writing: %x\n", RLE_indicator);

					// Write the length of 0 length region: Note that i does not have a 0, it must not be included.
					unsigned int l_RLE = (i - cur_0_run_start);
					fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
					fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
				} // RLE dump check.
				else
				{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
					fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

					// Save the position with normal dump where we dump 0's.
					for(int j = cur_0_run_start; j < i; j++)
					{
						// Write 0 to each position.
						unsigned char RLE_indicator = 0;
						fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
					} // j loop.
				} // non-RLE check.

				// Reset the RLE start position.
				cur_0_run_start = -1;
			} // 0-run check.	

			if(max_val == 0)
			{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
				fprintf(stderr, "Sanity check failed: Max value is 0 on a non-zero dumping path.\n");
				exit(0);
			}

			// Continue dumping the value for the current position.
			// Set the data range flag.
			if(max_val >= (1 << 8)-1)
			{
				drange_flag = 1;
			}

			// If val is not char, flag it to make sure we dump the right number of bytes.
			existence_flag = existence_flag | (drange_flag << 5);

			if(existence_flag == 0xFF)
			{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
				fprintf(stderr, "Sanity check failed: Existence flag is 0xFF.\n");
				exit(0);
			}

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Normal dumping value %d @ %d, Existence flag: 0x%x, drange_flag: %d\n", total_val, i, existence_flag, drange_flag);

			// Do a normal dump.
			// Write existence flag.
			fwrite(&existence_flag, sizeof(unsigned char), 1, f_op);

			// Dump the values.
			for(int part_i = 0; part_i < n_partitions; part_i++)
			{
				for(int allele_i = 0; allele_i < 5; allele_i++)
				{
					// Following ensures that we write only the positions that exist.
					//if(pileup[part_i][allele_i][i] > 0)
					if((existence_flag & (1 << allele_i)) > 0)
					{
						if(drange_flag == 1)
						{	
							unsigned short current_count_short = (unsigned short)(pileup[part_i][allele_i][i]);
							fwrite(&(current_count_short), sizeof(unsigned short), 1, f_op);
						}
						else
						{
							unsigned char current_count_char = (unsigned char)(pileup[part_i][allele_i][i]);
							fwrite(&(current_count_char), sizeof(unsigned char), 1, f_op);
						}
						//fprintf(stderr, "%d, ", pileup[part_i][allele_i][i]);
					} // count check for current position.
				} // allele_i loop.
			} // part_i loop.

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "\n");
		} // total value check.		

		// Get some simple stats.
		if(max_val == 0)
		{
			n_0_signal++;
		}
		else if(max_val <= ((2<<4)-1))
		{
			n_4bit_posns++;
		}
		else if(max_val <= ((2<<8)-1))
		{
			n_8bit_posns++;
		}
		else if(max_val <= ((2<<12)-1))
		{
			n_12bit_posns++;
		}
		else if(max_val <= ((2<<14)-1))
		{
			n_14bit_posns++;
		}
		else if(max_val <= ((2<<16)-1))
		{
			n_16bit_posns++;
		}
	} // i loop.

	// If there was a 0-run at the end of the file, dump the final RLE.
	if(cur_0_run_start != -1)
	{
		int i = l_sig+1;
		bool do_RLE_dump = ((i - cur_0_run_start) > 5);

		if(do_RLE_dump)
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

			// Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
			// Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
			unsigned char RLE_indicator = 0xFF;
			fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Writing: %x\n", RLE_indicator);

			// Write the length of 0 length region: Note that i does not have a 0, it must not be included.
			unsigned int l_RLE = (i - cur_0_run_start);
			fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
		} // RLE dump check.
		else
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

			// Save the position with normal dump where we dump 0's.
			for(int j = cur_0_run_start; j < i; j++)
			{
				// Write 0 to each position.
				unsigned char RLE_indicator = 0;
				fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
			} // j loop.
		} // non-RLE check.
	} // RLE check for the end of file.

	close_f(f_op, op_fp);
	fprintf(stderr, "%d 0 signal, %d 4-bit, %d 8-bit, %d 12-bit, %d 14-bit, %d 16-bit positions.\n", n_0_signal, n_4bit_posns, n_8bit_posns, n_12bit_posns, n_14bit_posns, n_16bit_posns);

	// Now, load and test if we received the correct values.
	
	bool loading_check = false;
	if(loading_check)
	{
		fprintf(stderr, "Loading and comparing.\n");

		// Load the compressed pileup; then do sanity check.
		int l_loaded_pileup = 0;
		int n_loaded_partitions = 0;
		unsigned short*** loaded_phred_partitioned_pileup = load_partitioned_compressed_pileups(op_fp, n_loaded_partitions, l_loaded_pileup);

		// Compare the loaded and the actual pileups.
		for(int part_i = 0; part_i < n_loaded_partitions; part_i++)
		{
			for(int i_allele = 0; i_allele < 5; i_allele++)
			{
				for(int i = 1; i <= l_sig; i++)
				{
					if(loaded_phred_partitioned_pileup[part_i][i_allele][i] != pileup[part_i][i_allele][i])
					{
						fprintf(stderr, "Sanity check failed: partition: %d; posn: %d; allele: %d: Loaded: %d, Original: %d\n", 
							part_i, i, i_allele, loaded_phred_partitioned_pileup[part_i][i_allele][i], pileup[part_i][i_allele][i]);
						exit(0);
					}
				} // i loop.				
			} // i_allele loop.
		} // part_i loop.

		// Free memory for both loaded and dumped pileups.
		for(int part_i = 0; part_i < n_partitions; part_i++)
		{
			for(int i_allele = 0; i_allele < 5; i_allele++)
			{
				delete [] loaded_phred_partitioned_pileup[part_i][i_allele];
			} // i_allele loop.
			delete [] loaded_phred_partitioned_pileup[part_i];
		} // part_i loop.
	}
	else
	{
		fprintf(stderr, "Skipping loading check.\n");
	}
}

void compress_nucleotide_pileup_track(unsigned short** pileup, int l_sig, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	
	// Write the length of pileup
	fwrite(&l_sig, sizeof(int), 1, f_op);
	int n_4bit_posns = 0;
	int n_8bit_posns = 0;
	int n_12bit_posns = 0;
	int n_14bit_posns = 0;
	int n_16bit_posns = 0;
	int n_0_signal = 0;

	// This is the 0 RLE start position.
	int cur_0_run_start = -1;
	for(int i = 1; i <= l_sig; i++)
	{
		// Dump the current 5 levels: 		
		// Generate the existence flag.
		//fprintf(stderr, "Position %d:\n", i);
		unsigned char existence_flag = 0;
		unsigned char drange_flag = 0;
		unsigned short max_val = 0;
		unsigned int total_val = 0;
		for(int allele_i = 0; allele_i < 5; allele_i++)
		{
			if(pileup[allele_i][i] > 0)
			{
				total_val += (unsigned int)(pileup[allele_i][i]);

				if(max_val < pileup[allele_i][i])
				{
					max_val = pileup[allele_i][i];
				}

				existence_flag = existence_flag | (1 << allele_i);
				//fprintf(stderr, "Nuc %d: %d\n", allele_i, pileup[allele_i][i]);
				//getc(stdin);
			} // signal check.
		} // allele_i loop.

		// Write the 1 byte phred-based partitioning flag.
		// Following is an RLE check.
		if(total_val == 0)
		{
			if(cur_0_run_start == -1)
			{
				// Initiate the 0 run.
				cur_0_run_start = i;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
				fprintf(stderr, "Initiated a 0-run @ %d\n", i);
			}
			else
			{
				// We are already in a zero run.
			}
		} // total_val 0 check.
		else
		{
			// Check if we will do a RLE dump.
			// Make sure there was a run of zeros before this position.
			if(cur_0_run_start != -1)
			{
				// Check the RLE length and dump if requested.
				bool do_RLE_dump = ((i - cur_0_run_start) > 5);

				if(do_RLE_dump)
				{
					//fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

					// Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
					// Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
					unsigned char RLE_indicator = 0xFF;
					fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
					fprintf(stderr, "Writing: %x\n", RLE_indicator);

					// Write the length of 0 length region: Note that i does not have a 0, it must not be included.
					unsigned int l_RLE = (i - cur_0_run_start);
					fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
					fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
				} // RLE dump check.
				else
				{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
					fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

					// Save the position with normal dump where we dump 0's.
					for(int j = cur_0_run_start; j < i; j++)
					{
						// Write 0 to each position.
						unsigned char RLE_indicator = 0;
						fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
					} // j loop.
				} // non-RLE check.

				// Reset the RLE start position.
				cur_0_run_start = -1;
			} // 0-run check.	

			if(max_val == 0)
			{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
				fprintf(stderr, "Sanity check failed: Max value is 0 on a non-zero dumping path.\n");
				exit(0);
			}
			
			// Set the data range flag: 1 for short, 0 for char.
			if(max_val >= (1 << 8)-1)
			{
				drange_flag = 1;
			}

			// If val is not char, flag it to make sure we dump the right number of bytes.
			existence_flag = existence_flag | (drange_flag << 5);

			// Write existence flag.
			fwrite(&existence_flag, sizeof(unsigned char), 1, f_op);

			// Dump the numbers.
			for(int allele_i = 0; allele_i < 5; allele_i++)
			{
				if(pileup[allele_i][i] > 0)
				{
					if(drange_flag == 1)
					{	
						unsigned short current_count_short = (unsigned short)(pileup[allele_i][i]);
						fwrite(&(current_count_short), sizeof(unsigned short), 1, f_op);
					}
					else
					{
						unsigned char current_count_char = (unsigned char)(pileup[allele_i][i]);
						fwrite(&(current_count_char), sizeof(unsigned char), 1, f_op);
					}
				} // count check for current position.
			} // allele_i loop.

		} // total value check.

		// Get some simple stats.
		if(max_val == 0)
		{
			n_0_signal++;
		}
		else if(max_val <= ((2<<4)-1))
		{
			n_4bit_posns++;
		}
		else if(max_val <= ((2<<8)-1))
		{
			n_8bit_posns++;
		}
		else if(max_val <= ((2<<12)-1))
		{
			n_12bit_posns++;
		}
		else if(max_val <= ((2<<14)-1))
		{
			n_14bit_posns++;
		}
		else if(max_val <= ((2<<16)-1))
		{
			n_16bit_posns++;
		}
	} // i loop.
	close_f(f_op, op_fp);
	fprintf(stderr, "%d 0 signal, %d 4-bit, %d 8-bit, %d 12-bit, %d 14-bit, %d 16-bit positions.\n", n_0_signal, n_4bit_posns, n_8bit_posns, n_12bit_posns, n_14bit_posns, n_16bit_posns);

	// Now, load and test if we received the correct values.
	bool loading_check = false;
	if(loading_check)
	{
		fprintf(stderr, "Loading and comparing.\n");

		// Load the compressed pileup; then do sanity check.
		int l_loaded_pileup = 0;
		unsigned short** loaded_pileup = load_compressed_pileups(op_fp, l_loaded_pileup);

		// Compare the loaded and the actual pileups.
		for(int i_allele = 0; i_allele < 5; i_allele++)
		{
			for(int i = 1; i <= l_sig; i++)
			{
				if(loaded_pileup[i_allele][i] != pileup[i_allele][i])
				{
					fprintf(stderr, "Sanity check failed: posn: %d; allele: %d: %d,%d\n", 
						i, i_allele, loaded_pileup[i_allele][i], pileup[i_allele][i]);
					exit(0);
				}
			} // i loop.				
		} // i_allele loop.

		// Free memory for both loaded and dumped pileups.
		for(int i_allele = 0; i_allele < 5; i_allele++)
		{
			delete [] loaded_pileup[i_allele];
		} // i_allele loop.
	}
	else
	{
		fprintf(stderr, "Skipping loading check.\n");
	}
}



void dump_nucleotide_pileup_per_SAM_file_orig(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, int& n_processed_reads)
{
	//// Init the file.
	//char summary_fp[1000];
	//sprintf(summary_fp, "%s/%s", op_dir, t_config_params::OP_filenames[OP_PILEUP_SUMMARY_FN]);
	//FILE* f_summary = open_f(summary_fp, "w");
	//fclose(f_summary);

	// Initialize the number of processed readss.
	n_processed_reads = 0;

	// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
	bool allele_counts_are_there = true;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
		if(!check_file(cur_allele_fp))
		{
			// This file does not exist, not all the allele counts are there.
			allele_counts_are_there = false;
			break;
		}
	} // i_chr loop.

	if(allele_counts_are_there)
	{
		fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
		return;
	}

	fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

	fprintf(stderr, "Allocating pileup memory.\n");
	unsigned short*** per_chrom_nuc_count_per_allele = new unsigned short**[(int)chr_ids->size()];
	//int** per_chrom_coverage = new int*[(int)chr_ids->size()];
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_sig = chr_lengths->at(i_chr);

		//per_chrom_coverage[i_chr] = new int[l_sig + 2];
		//memset(per_chrom_coverage[i_chr], 0, sizeof(int) * l_sig);

		per_chrom_nuc_count_per_allele[i_chr] = allocate_pileup(l_sig);
	} // i_chr loop.
	fprintf(stderr, "Done.\n");

	int max_n_alleles_per_posn = 1<<16;
	fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

	// Start reading the file.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = fopen(sam_fp, "r");
	}

	int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
	//for(int i = 0; i < 256; i++)
	//{
	//	per_nuc_val[i] = 4;
	//} // i loop.
	//per_nuc_val[(int)'A'] = 0;
	//per_nuc_val[(int)'C'] = 1;
	//per_nuc_val[(int)'G'] = 2;
	//per_nuc_val[(int)'T'] = 3;
	//per_nuc_val[(int)'U'] = 3;
	//per_nuc_val[(int)'a'] = 0;
	//per_nuc_val[(int)'c'] = 1;
	//per_nuc_val[(int)'g'] = 2;
	//per_nuc_val[(int)'t'] = 3;
	//per_nuc_val[(int)'u'] = 3;

	// Enter file reading loop.
	char read_id[1000];
	//int flag;
	char chrom[100];
	int chr_index;
	char mapping_map_str[10000];
	char cur_fragment[1000];
	char flag_str[100];
	char _chr_index_str[100];
	char phred_quality_str[100000];
	char mapp_quality_str[100];

	int n_unmapped_reads = 0;
	int n_low_quality_reads = 0;

	int phred_score_base = -123;

	while(1)
	{
		char* cur_line = getline(f_sam);

		if(cur_line == NULL)
		{
			break;
		}

		// If this is a comment line, skip it.
		if(cur_line[0] == '@')
		{
			delete [] cur_line;
			continue;
		}

		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
		{
			//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
			// If the quality is not adequate, do not use this read.
			if(atoi(mapp_quality_str) < min_mapp_qual)
			{
				n_low_quality_reads++;
				delete [] cur_line;
				continue;
			}

			// Make sure the normalized chromosome ids match.
			normalize_chr_id(chrom);

			int flag = atoi(flag_str);

			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If we do not have the chromosome, do not process.
			if(i_chr == (int)chr_ids->size())
			{
				n_unmapped_reads++;
				delete [] cur_line;
				continue;
			}

			int _chr_index = atoi(_chr_index_str);

			// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
			chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

			// Sanity check. Is this fragment mapped?
			if(flag & 0x04)
			{
				// The read is not mapping.
			}
			else
			{
				n_processed_reads++;
					
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %d. read             \r", n_processed_reads);
				}

				int i_mapp_map = 0;
				//t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
				while(mapping_map_str_valid && 
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map, 
														is_matching,
														l_cur_entry,
														entry_type_char);

					// If this block is matching, update the pileup.
					if(is_matching)
					{
						// Update the counts for the nucleotides.
						int cur_read_nuc_i = read_nuc_index;
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							//int cur_nuc_num = nuc_2_num(cur_fragment[cur_read_nuc_i]);
							int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

							//if(numerized_sequence_signal[cur_genome_i] > 0 &&
							//	nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
							//{
							//	fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
							//}

							if(cur_nuc_num < 4)
							{
								// Update the count: The allelic counts must be checked for bounds.
								if(per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
								{
									per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i]++;
								}

								//per_chrom_coverage[i_chr][cur_genome_i]++;
								cur_read_nuc_i++;

								//fprintf(stderr, "%d: %c, %c\n", i_cur, num_2_nuc(cur_nuc_num), num_2_nuc(numerized_sequence_signal[i_cur]));
							}
						} // i_cur loop.
					}
					else if(entry_type_char == 'D')
					{
						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							per_chrom_nuc_count_per_allele[i_chr][4][cur_genome_i]++;
						} // cur_genome_i loop.
					}
					else if(entry_type_char == 'I')
					{
						// Insertion to the reference: This is included to one position.
						per_chrom_nuc_count_per_allele[i_chr][4][chr_index]++;
					}
					// Update the base for the current entry.
					if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				//delete(cur_entry_length_str);
			} // mapping check for the current read.
		} // SAM read line parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "Finished reading.\n");

	fclose(f_sam);

	// Dump the counts per position.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));

		// compress and dump the current file.
		compress_nucleotide_pileup_track(per_chrom_nuc_count_per_allele[i_chr], chr_lengths->at(i_chr), cur_allele_fp);
	} // i_chr loop.

	//char analysis_summary_str[1000];
	//sprintf(analysis_summary_str, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality, n_unmapped_reads);
	//t_config_params::copy_analysis_summary_string(analysis_summary_str);

	//// Write the summary file.
	//f_summary = open_f(summary_fp, "w");
	//fprintf(f_summary, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality_reads, n_unmapped_reads);
	//fclose(f_summary);
}

void dump_nucleotide_pileup_per_SAM_file(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, unsigned long long& n_processed_reads)
{
	fprintf(stderr, "Generating pileup from SAM file with min mapp qual %d, min base qual %d\n", min_mapp_qual, min_phred_qual);

	//// Init the file.
	//char summary_fp[1000];
	//sprintf(summary_fp, "%s/%s", op_dir, t_config_params::OP_filenames[OP_PILEUP_SUMMARY_FN]);
	//FILE* f_summary = open_f(summary_fp, "w");
	//fclose(f_summary);

	// Initialize the number of processed readss.
	n_processed_reads = 0;

	//// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
	//bool allele_counts_are_there = true;
	//for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	//{
	//	char cur_allele_fp[1000];
	//	sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
	//	if(!check_file(cur_allele_fp))
	//	{
	//		// This file does not exist, not all the allele counts are there.
	//		allele_counts_are_there = false;
	//		break;
	//	}
	//} // i_chr loop.

	//if(allele_counts_are_there)
	//{
	//	fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
	//	return;
	//}

	fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

	fprintf(stderr, "Allocating pileup memory.\n");
	unsigned short*** per_chrom_nuc_count_per_allele = new unsigned short**[(int)chr_ids->size()];
	//int** per_chrom_coverage = new int*[(int)chr_ids->size()];
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_sig = chr_lengths->at(i_chr);

		// Check the existing pileup file, if it is there, add to it.
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
		if(check_file(cur_allele_fp))
		{
			int l_loaded_pileup = 0;
			per_chrom_nuc_count_per_allele[i_chr] = load_compressed_pileups(cur_allele_fp, l_loaded_pileup);
			if(l_loaded_pileup != chr_lengths->at(i_chr))
			{
				fprintf(stderr, "Loaded pileup is not the same length as the chromosome info: %s: %d, %d\n", chr_ids->at(i_chr), chr_lengths->at(i_chr), l_loaded_pileup);
				exit(0);
			}
			else
			{
				fprintf(stderr, "Loaded existing pileup from %s.\n", cur_allele_fp);
			}
		}
		else
		{
			per_chrom_nuc_count_per_allele[i_chr] = allocate_pileup(l_sig);
		}
	} // i_chr loop.
	fprintf(stderr, "Done.\n");

	int max_n_alleles_per_posn = 1<<16;
	fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

	// Start reading the file.
	FILE* f_sam = open_f(sam_fp, "r");

	int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
	//for(int i = 0; i < 256; i++)
	//{
	//	per_nuc_val[i] = 4;
	//} // i loop.
	//per_nuc_val[(int)'A'] = 0;
	//per_nuc_val[(int)'C'] = 1;
	//per_nuc_val[(int)'G'] = 2;
	//per_nuc_val[(int)'T'] = 3;
	//per_nuc_val[(int)'U'] = 3;
	//per_nuc_val[(int)'a'] = 0;
	//per_nuc_val[(int)'c'] = 1;
	//per_nuc_val[(int)'g'] = 2;
	//per_nuc_val[(int)'t'] = 3;
	//per_nuc_val[(int)'u'] = 3;

	// Enter file reading loop.
	char read_id[1000];
	//int flag;
	char chrom[100];
	int chr_index;
	char mapping_map_str[10000];
	char cur_fragment[1000];
	char flag_str[100];
	char _chr_index_str[100];
	char phred_quality_str[100000];
	char mapp_quality_str[100];

	int n_unmapped_reads = 0;
	int n_low_quality_reads = 0;

	int phred_score_base = -123;

	while(1)
	{
		char* cur_line = getline(f_sam);

		if(cur_line == NULL)
		{
			break;
		}

		// If this is a comment line, skip it.
		if(cur_line[0] == '@')
		{
			delete [] cur_line;
			continue;
		}

		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
		{
			bool phred_exists = true;
			if (t_string::compare_strings(phred_quality_str, "*") ||
				t_string::string_length(phred_quality_str) != t_string::string_length(cur_fragment))
			{
				phred_exists = false;
			}

			// Note that this can be updated for any read that contains a valid quality score.
			if(phred_exists &&
				phred_score_base == -123)
			{
				int l_read = strlen(phred_quality_str);
				for(int i = 0; i < l_read; i++)
				{
					if(phred_quality_str[i] > 'J')
					{
						fprintf(stderr, "Phred+64 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = ';';
						break;
					}
					else if(phred_quality_str[i] < ';')
					{
						fprintf(stderr, "Phred+33 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = '!';
						break;
					}
				} // i loop.
			}
			
			//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
			// If the quality is not adequate, do not use this read.
			if(atoi(mapp_quality_str) < min_mapp_qual)
			{
				n_low_quality_reads++;
				delete [] cur_line;
				continue;
			}

			// Make sure the normalized chromosome ids match.
			normalize_chr_id(chrom);

			int flag = atoi(flag_str);

			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If we do not have the chromosome, do not process.
			if(i_chr == (int)chr_ids->size())
			{
				n_unmapped_reads++;
				delete [] cur_line;
				continue;
			}

			int _chr_index = atoi(_chr_index_str);

			// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
			chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

			// Sanity check. Is this fragment mapped?
			if(flag & 0x04)
			{
				// The read is not mapping.
			}
			else
			{
				n_processed_reads++;
					
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %llu. read             \r", n_processed_reads);
				}

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				fprintf(stderr, "Processing:\n%s\n%s\n", cur_fragment, phred_quality_str);
}

				int i_mapp_map = 0;
				//t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
				while(mapping_map_str_valid && 
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map, 
														is_matching,
														l_cur_entry,
														entry_type_char);

					// If this block is matching, update the pileup.
					if(is_matching)
					{
						// Update the counts for the nucleotides.
						int cur_read_nuc_i = read_nuc_index;
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

							//if(numerized_sequence_signal[cur_genome_i] > 0 &&
							//	nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
							//{
							//	fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
							//}

							if(cur_nuc_num < 4)
							{
								// Update the count: The allelic counts must be checked for bounds.
								if(per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
								{
									// If phred is missing, process all of them, if not, check the quality.
									bool phred_check = true;
									if (phred_exists)
									{
										phred_check = ((phred_quality_str[cur_read_nuc_i] - phred_score_base) > min_phred_qual);
									}

									// Does phred check hold?
									if(phred_check)
									{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
										fprintf(stderr, "Adding %d: %c, %c (%d)\n", 
												cur_read_nuc_i, 
												cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
												(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}

										// Phred quality holds, update the pileup position.
										if (cur_genome_i <= chr_lengths->at(i_chr))
										{
											per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i]++;
										}
										else
										{
											fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
												chr_ids->at(i_chr), cur_genome_i, chr_lengths->at(i_chr));
										}
									} // phred qual pass check.
									else
									{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
										fprintf(stderr, "Skipping %d: %c, %c (%d)\n", 
											cur_read_nuc_i, 
											cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
											(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}
									} // phred qual nopass check.
								} // max_n_alleles_per_posn check.
							} // char < 4 check.

							// Update the read's nucleotide index.
							cur_read_nuc_i++;
						} // i_cur loop.
					} // genomic match check.
					else if(entry_type_char == 'D')
					{
						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							if (cur_genome_i <= chr_lengths->at(i_chr))
							{
								per_chrom_nuc_count_per_allele[i_chr][4][cur_genome_i]++;
							}
							else
							{
								fprintf(stderr, "%s:%d is further away from what we have in the lengths of the chromosomes (%d); can there be a mismatch between assembly that reads are mapped to vs this?\n",
									chr_ids->at(i_chr), cur_genome_i, chr_lengths->at(i_chr));
							}
						} // cur_genome_i loop.
					}
					else if(entry_type_char == 'I')
					{
						// Insertion to the reference: This is included to one position.
						if (chr_index < chr_lengths->at(i_chr))
						{
							per_chrom_nuc_count_per_allele[i_chr][4][chr_index]++;
						}
					}
					// Update the base for the current entry.
					if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				//delete(cur_entry_length_str);
			} // mapping check for the current read.
		} // SAM read line parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "Finished reading.\n");

	close_f(f_sam, sam_fp);

	// Dump the counts per position.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin.gz", op_dir, chr_ids->at(i_chr));

		// compress and dump the current file.
		compress_nucleotide_pileup_track(per_chrom_nuc_count_per_allele[i_chr], chr_lengths->at(i_chr), cur_allele_fp);
	} // i_chr loop.

	//char analysis_summary_str[1000];
	//sprintf(analysis_summary_str, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality, n_unmapped_reads);
	//t_config_params::copy_analysis_summary_string(analysis_summary_str);

	//// Write the summary file.
	//f_summary = open_f(summary_fp, "w");
	//fprintf(f_summary, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality_reads, n_unmapped_reads);
	//fclose(f_summary);
}

void dump_nucleotide_pileup_per_SAM_file_phred_partitioning(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, vector<int>* phred_qual_partitions, unsigned long long& n_processed_reads)
{
	// Sort the quality partition boundary values.
	int max_qual = 10000;
	phred_qual_partitions->push_back(0);
	phred_qual_partitions->push_back(max_qual);
	sort(phred_qual_partitions->begin(), phred_qual_partitions->end());
	fprintf(stderr, "Generating pileup from SAM file with min mapp qual %d, with %d phred partitions @:\n", min_mapp_qual, phred_qual_partitions->size()-1);

	// Check the quality partitions and merge the ones that overlap.
	vector<int>* uniq_phred_qual_partitions = new vector<int>();
	for(int i_p = 0; i_p < phred_qual_partitions->size(); i_p++)
	{
		if(i_p+1 < phred_qual_partitions->size() &&
			phred_qual_partitions->at(i_p) == phred_qual_partitions->at(i_p+1))
		{
			fprintf(stderr, "Merging partitions that match at qualities of %d\n", phred_qual_partitions->at(i_p));
		}
		else
		{
			uniq_phred_qual_partitions->push_back(phred_qual_partitions->at(i_p));
		}
	} // i_p loop.

	// Replace the phred quality partitions.
	phred_qual_partitions->clear();
	phred_qual_partitions->insert(phred_qual_partitions->end(), uniq_phred_qual_partitions->begin(), uniq_phred_qual_partitions->end());

	int* part_i_per_phred_qual = new int[max_qual+1];
	for(int i_p = 1; i_p < phred_qual_partitions->size(); i_p++)
	{
		fprintf(stderr, "Partition %d: [%d-%d]\n", i_p-1, phred_qual_partitions->at(i_p-1)+1, phred_qual_partitions->at(i_p));
		for(int i = phred_qual_partitions->at(i_p-1)+1; i <= phred_qual_partitions->at(i_p); i++)
		{
			part_i_per_phred_qual[i] = i_p - 1;
		} // i loop.
	} // i_p loop.
	part_i_per_phred_qual[0] = 0;

	int n_partitions = phred_qual_partitions->size()-1;

	// Initialize the number of processed readss.
	n_processed_reads = 0;

	// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
	bool allele_counts_are_there = true;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
		if(!check_file(cur_allele_fp))
		{
			// This file does not exist, not all the allele counts are there.
			allele_counts_are_there = false;
			break;
		}
	} // i_chr loop.

	if(allele_counts_are_there)
	{
		fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
		return;
	}

	fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

	fprintf(stderr, "Allocating pileup memory.\n");
	unsigned short**** per_chrom_per_part_nuc_count_per_allele = new unsigned short***[(int)chr_ids->size()];
	double total_mem = 0;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		int l_sig = chr_lengths->at(i_chr);

		per_chrom_per_part_nuc_count_per_allele[i_chr] = new unsigned short**[n_partitions];

		for(int part_i = 0; part_i < n_partitions; part_i++)
		{
			per_chrom_per_part_nuc_count_per_allele[i_chr][part_i] = allocate_pileup(l_sig);
		} // part_i loop.
	} // i_chr loop.
	fprintf(stderr, "\nDone.\n");

	int max_n_alleles_per_posn = 1<<16;
	fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

	// Start reading the file.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = fopen(sam_fp, "r");
	}

	int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
	//for(int i = 0; i < 256; i++)
	//{
	//	per_nuc_val[i] = 4;
	//} // i loop.
	//per_nuc_val[(int)'A'] = 0;
	//per_nuc_val[(int)'C'] = 1;
	//per_nuc_val[(int)'G'] = 2;
	//per_nuc_val[(int)'T'] = 3;
	//per_nuc_val[(int)'U'] = 3;
	//per_nuc_val[(int)'a'] = 0;
	//per_nuc_val[(int)'c'] = 1;
	//per_nuc_val[(int)'g'] = 2;
	//per_nuc_val[(int)'t'] = 3;
	//per_nuc_val[(int)'u'] = 3;

	// Enter file reading loop.
	char read_id[1000];
	//int flag;
	char chrom[100];
	int chr_index;
	char mapping_map_str[10000];
	char cur_fragment[1000];
	char flag_str[100];
	char _chr_index_str[100];
	char phred_quality_str[100000];
	char mapp_quality_str[100];

	int n_unmapped_reads = 0;
	int n_low_quality_reads = 0;

	int phred_score_base = -123;

	while(1)
	{
		char* cur_line = getline(f_sam);

		if(cur_line == NULL)
		{
			break;
		}

		// If this is a comment line, skip it.
		if(cur_line[0] == '@')
		{
			delete [] cur_line;
			continue;
		}

		if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
		{
			if(phred_score_base == -123)
			{
				int l_read = strlen(phred_quality_str);
				for(int i = 0; i < l_read; i++)
				{
					if(phred_quality_str[i] > 'J')
					{
						fprintf(stderr, "Phred+64 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = ';';
						break;
					}
					else if(phred_quality_str[i] < ';')
					{
						fprintf(stderr, "Phred+33 encoding @ %llu. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
							cur_fragment, phred_quality_str);
						phred_score_base = '!';
						break;
					}
				} // i loop.
			}

			//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
			// If the quality is not adequate, do not use this read.
			if(atoi(mapp_quality_str) < min_mapp_qual)
			{
				n_low_quality_reads++;
				delete [] cur_line;
				continue;
			}

			// Make sure the normalized chromosome ids match.
			normalize_chr_id(chrom);

			int flag = atoi(flag_str);

			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If we do not have the chromosome, do not process.
			if(i_chr == (int)chr_ids->size())
			{
				n_unmapped_reads++;
				delete [] cur_line;
				continue;
			}

			int _chr_index = atoi(_chr_index_str);

			// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
			chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

			// Sanity check. Is this fragment mapped?
			if(flag & 0x04)
			{
				// The read is not mapping.
			}
			else
			{
				n_processed_reads++;
					
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %llu. read             \r", n_processed_reads);
				}

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				fprintf(stderr, "Processing:\n%s\n%s\n", cur_fragment, phred_quality_str);
}

				int i_mapp_map = 0;
				//t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				int read_nuc_index = 0;

				// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
				while(mapping_map_str_valid && 
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
														i_mapp_map, 
														is_matching,
														l_cur_entry,
														entry_type_char);

					// If this block is matching, update the pileup.
					if(is_matching)
					{
						// Update the counts for the nucleotides.
						int cur_read_nuc_i = read_nuc_index;
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

							//if(numerized_sequence_signal[cur_genome_i] > 0 &&
							//	nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
							//{
							//	fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
							//}

							if(cur_nuc_num < 4)
							{
								int cur_nuc_phred_qual = phred_quality_str[cur_read_nuc_i] - phred_score_base;

								// Get the partition index.
								int cur_read_phred_part_i = part_i_per_phred_qual[cur_nuc_phred_qual];

								// Update the count: The allelic counts must be checked for bounds.
								if(per_chrom_per_part_nuc_count_per_allele[i_chr][cur_read_phred_part_i][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
								{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
									fprintf(stderr, "Adding %d: %c, %c (%d)\n", 
											cur_read_nuc_i, 
											cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
											(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}

										// Phred quality holds, update the pileup position.
										per_chrom_per_part_nuc_count_per_allele[i_chr][cur_read_phred_part_i][cur_nuc_num][cur_genome_i]++;
//									} // phred qual pass check.
//									else
//									{
//if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
//{
//										fprintf(stderr, "Skipping %d: %c, %c (%d)\n", 
//											cur_read_nuc_i, 
//											cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
//											(int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
//}
//									} // phred qual nopass check.
								} // max_n_alleles_per_posn check.
							} // char < 4 check.

							// Update the read's nucleotide index.
							cur_read_nuc_i++;
						} // i_cur loop.
					} // genomic match check.
					else if(entry_type_char == 'D')
					{
						// Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
						for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
						{
							per_chrom_per_part_nuc_count_per_allele[i_chr][n_partitions-1][4][cur_genome_i]++;
						} // cur_genome_i loop.
					}
					else if(entry_type_char == 'I')
					{
						// Insertion to the reference: This is included to one position.
						per_chrom_per_part_nuc_count_per_allele[i_chr][n_partitions-1][4][chr_index]++;
					}
					// Update the base for the current entry.
					if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					// Update the base for the current read if requested.
					if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					{
						read_nuc_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				//delete(cur_entry_length_str);
			} // mapping check for the current read.
		} // SAM read line parse check.
		else
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "Finished reading.\n");

	fclose(f_sam);

	// Dump the counts per position.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_allele_fp[1000];
		sprintf(cur_allele_fp, "%s/%s_allele_counts.bin.gz", op_dir, chr_ids->at(i_chr));

		// compress and dump the current file.
		compress_partitioned_nucleotide_pileup_track(per_chrom_per_part_nuc_count_per_allele[i_chr], n_partitions, chr_lengths->at(i_chr), cur_allele_fp);
	} // i_chr loop.
}

void preprocess_SAM_read_line_positional_info(char* cur_line,
	char* chrom,
	int& chr_index, 
	int& mapQ,
	char& strand_char,
	char* cigar_str)
{
	// Skip the comment and headers.
	if (cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	char flag_str[100];
	int _chr_index;
	char _chr_index_str[100];
	char mapQ_str[1000];

	int char_i = 0;
	if (!t_string::get_next_token(cur_line, NULL, 0, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	if (!t_string::get_next_token(cur_line, flag_str, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	if (!t_string::get_next_token(cur_line, chrom, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	if (!t_string::get_next_token(cur_line, _chr_index_str, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}
	
	if (!t_string::get_next_token(cur_line, mapQ_str, 1000, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}	

	if(!t_string::get_next_token(cur_line, cigar_str, 100, "\t", char_i))
	{
		chrom[0] = 0;
		return;
	}

	// Is this fragment mapped?
	if (flag & 0x04)
	{
		// The read is not mapping.
		chrom[0] = 0;
		mapQ = 0;
		chr_index = 0;
		flag = 0;
		strand_char = 0;
	}
	else
	{
		// Tramslate the coordinate and strand.
		_chr_index = atoi(_chr_index_str);
		mapQ = atoi(mapQ_str);
		chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);
		flag = atoi(flag_str);

		strand_char = 'F';
		if (flag & 0x10)
		{
			strand_char = 'R';
		}
	}

	//t_string::clean_tokens(cur_tokens);
}


void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	char flag_str[100];
	int _chr_index;
	char _chr_index_str[100];
	char fragment[100000];
	char phred_quality_str[100000];

	//t_string_tokens* cur_tokens = t_string::tokenize_by_chars(cur_line, "\t");
	//if(sscanf(cur_line, "%s %d %s %d %*s %s %*s %*s %*s %s %s", read_id, &flag, chrom, &_chr_index, cigar_str, fragment, phred_quality_str) == 7)
	//if(cur_tokens->size() >= 11)
	if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, cigar_str, fragment, phred_quality_str) == 7)
	{
		//t_string::copy(read_id, cur_tokens->at(0)->str());
		//flag = atoi(cur_tokens->at(1)->str());
		//t_string::copy(chrom, cur_tokens->at(2)->str());
		//_chr_index = atoi(cur_tokens->at(3)->str());
		//t_string::copy(cigar_str, cur_tokens->at(5)->str());
		//t_string::copy(fragment, cur_tokens->at(9)->str());
		//t_string::copy(phred_quality_str, cur_tokens->at(10)->str());

		_chr_index = atoi(_chr_index_str);
		flag = atoi(flag_str);

		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(flag & 0x10)
		{
			strand_char = 'R';
		}

		// Sanity check. Is this fragment mapped?
		if(flag & 0x04)
		{
			// The read is not mapping.
			chrom[0] = 0;
		}
		else
		{
			chr_index = _chr_index;
			sequenced_length = strlen(fragment);
		}
	}
	else
	{
		chrom[0] = 0;
	}

	//t_string::clean_tokens(cur_tokens);
}

void preprocess_PE_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	bool& first_segment_in_template,
	bool& last_segment_in_template,
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	int& mapping_quality,
	char* cigar_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	int _chr_index = 0;
	int _mapping_quality = 0;
	char fragment[100000];
	char phred_quality_str[100000];

	first_segment_in_template = false;
	last_segment_in_template = false;

	if(sscanf(cur_line, "%s %d %s %d %d %s %*s %*s %*s %s %s", read_id, &flag, chrom, &_chr_index, &_mapping_quality, cigar_str, fragment, phred_quality_str) == 8)
	{
		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(flag & 0x10)
		{
			strand_char = 'R';
		}

		// Sanity check. Is this fragment mapped?
		if(flag & 0x04)
		{
			// The read is not mapping.
			chrom[0] = 0;
		}
		else
		{
			if(flag & 0x40)
			{
				first_segment_in_template = true;
			}
			else if(flag & 0x80)
			{
				last_segment_in_template = true;
			}
			else if(flag & 0x01) // If the read has multiple fragments, this is a problem.
			{
				fprintf(stderr, "The segment is not first or last segment in the read that has multiple fragments: %s: %d\n", read_id, flag);
				exit(0);
			}

			chr_index = _chr_index;
			mapping_quality = _mapping_quality;
			sequenced_length = strlen(fragment);
		}
	}
	else
	{
		chrom[0] = 0;
	}
}


//void preprocess_BED_read_line(char* cur_line, 
//	char* chrom, 
//	int& chr_index, int& sequenced_length, 
//	char& strand_char, 
//	char* mapping_quality_str)
//{
//	// Skip the comment and headers.
//	if(cur_line[0] == '@')
//	{
//		chrom[0] = 0;
//		chr_index = 0;
//		return;
//	}
//
//	int _chr_start;
//	int _chr_end;
//	char strand;
//
//	if(sscanf(cur_line, "%s %d %d %*s %*s %c", chrom, &_chr_start, &_chr_end, &strand) == 6)
//	{
//		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
//		int chr_start = translate_coord(_chr_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
//		int chr_end = translate_coord(_chr_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
//
//		// Check the flag and determine the strand.
//		strand_char = 'F';
//		if(strand_char == '-')
//		{
//			strand_char = 'R';
//		}
//
//		chr_index = chr_start;
//		sequenced_length = chr_end-chr_start+1;
//	}
//	else
//	{
//		chrom[0] = 0;
//	}
//}

void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char cur_fragment[100];
	char quality_str[100];
	int _chr_index;
	char _strand_char;

	if(sscanf(cur_line, "%s %s %s %*d %*d %*d %s %d %c", read_id, cur_fragment, quality_str, chrom, &_chr_index, &_strand_char) == 6)
	{
		chr_index = _chr_index;
		sequenced_length = strlen(cur_fragment);
		sprintf(mapping_quality_str, "%dM", sequenced_length);

		strand_char = 'F';
		if(_strand_char == '-')
		{
			strand_char = 'R';
		}
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	char nucs[1000];
	if(sscanf(cur_line, "%s %c %s %d %s", read_id, &strand_sign, chrom, &chr_start_index, nucs) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		chr_start_index = translate_coord(chr_start_index, BOWTIE_COORDS::start_base, CODEBASE_COORDS::start_base);

		sprintf(mapping_quality_str, "%dM", (int)strlen(nucs));

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = strlen(nucs);
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;

	const int l_buff = 1000;
	char per_entry_buff[4][l_buff];
	int i_cur_char = 0;

	if (!t_string::get_next_token(cur_line, per_entry_buff[0], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[1], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[2], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[3], 1000, "\t", i_cur_char))
	{
		chrom[0] = 0;
		return;
	}
	else
	{
		strcpy(chrom, per_entry_buff[0]);
		chr_start_index = atoi(per_entry_buff[1]);
		chr_end_index = atoi(per_entry_buff[2]);
		strand_sign = per_entry_buff[3][0];
	}

	// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
	sprintf(mapping_quality_str, "%dM", chr_end_index - chr_start_index);
	sequenced_length = chr_end_index - chr_start_index;

	chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

	// Check the flag and determine the strand.
	strand_char = 'F';
	if (strand_sign == '-')
	{
		strand_char = 'R';
	}

	chr_index = chr_start_index;
}

void preprocess_BED5_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;

	const int l_buff = 1000;
	char per_entry_buff[5][l_buff];
	int i_cur_char = 0;

	if (!t_string::get_next_token(cur_line, per_entry_buff[0], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[1], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[2], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[3], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[4], 1000, "\t", i_cur_char))
	{
		chrom[0] = 0;
		return;
	}
	else
	{
		strcpy(chrom, per_entry_buff[0]);
		chr_start_index = atoi(per_entry_buff[1]);
		chr_end_index = atoi(per_entry_buff[2]);
		strand_sign = per_entry_buff[4][0];
	}

	// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
	sprintf(mapping_quality_str, "%dM", chr_end_index - chr_start_index);
	sequenced_length = chr_end_index - chr_start_index;

	chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

	// Check the flag and determine the strand.
	strand_char = 'F';
	if (strand_sign == '-')
	{
		strand_char = 'R';
	}

	chr_index = chr_start_index;
}

void preprocess_BED6_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;

	const int l_buff = 1000;
	char per_entry_buff[6][l_buff];
	int i_cur_char = 0;

	if (!t_string::get_next_token(cur_line, per_entry_buff[0], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[1], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[2], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[3], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[4], 1000, "\t", i_cur_char) ||
		!t_string::get_next_token(cur_line, per_entry_buff[5], 1000, "\t", i_cur_char))
	{
		chrom[0] = 0;
		return;
	}
	else
	{
		strcpy(chrom, per_entry_buff[0]);
		chr_start_index = atoi(per_entry_buff[1]);
		chr_end_index = atoi(per_entry_buff[2]);
		strand_sign = per_entry_buff[5][0];
	}

	// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
	sprintf(mapping_quality_str, "%dM", chr_end_index - chr_start_index);
	sequenced_length = chr_end_index - chr_start_index;

	chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

	// Check the flag and determine the strand.
	strand_char = 'F';
	if (strand_sign == '-')
	{
		strand_char = 'R';
	}

	chr_index = chr_start_index;
}

void preprocess_BED456_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str)
{
	int n_cols = t_string::fast_count_tokens(cur_line, true, "\t");
	
	/*fprintf(stderr, "%d columns.\r", n_cols);*/

	if (n_cols == 6)
	{
		preprocess_BED6_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);
	}
	else if(n_cols == 5)
	{
		preprocess_BED5_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);
	}
	else if (n_cols == 4)
	{
		preprocess_BED4_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			mapping_quality_str);
	}
}

double get_n_mapped_nucs(vector<t_mapped_fragment*>* fragments)
{
	double n_mapped_nucs = 0.0;

	for(int i_frag = 0; i_frag < (int)fragments->size(); i_frag++)
	{
		n_mapped_nucs += fragments->at(i_frag)->sequenced_fragment_length;
	} // i_frag loop.

	return(n_mapped_nucs);
}

void buffer_per_nucleotide_preprocessed_read_profile_no_buffer(char* sorted_read_fp,
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int max_l_read,
	int l_buffer, int& l_data)
{
	for(int i = 0; i <= l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initialize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	char strand_char = 0;
	//char cur_fragment[10000];
	char mapping_map_str[10000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = open_f(sorted_read_fp, "r");
	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
		}	

		n_processed_reads++;

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		int read_start = 0;
		int read_end = 0;

		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			if(is_matching)
			{
				if(read_start == 0)
				{
					read_start = chr_index;
				}
			} // CIGAR string matching check

			read_end = chr_index + l_cur_entry - 1;

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Increase the height: The indexing is already 1 based, there is no conversion needed.
		if(signal_profile_buffer != NULL &&
			max_l_read >= (read_end-read_start+1))
		{
			for(int i_cur = read_start; i_cur <= read_end; i_cur++)
			{
				signal_profile_buffer[i_cur]++;
			} // i_cur loop.
		}

		// Update the strand signals, if requested.
		if(forward_strand_signal != NULL &&
			max_l_read >= (read_end-read_start+1))
		{
			if(strand_char == 'F')
			{
				// Update the forward strand signal.
				for(int i_cur = read_start; i_cur <= read_end; i_cur++)
				{
					forward_strand_signal[i_cur]++;
				} // i_cur loop.
			}
			else
			{
				// Update the reverse strand signal.
				for(int i_cur = read_start; i_cur <= read_end; i_cur++)
				{
					reverse_strand_signal[i_cur]++;
				} // i_cur loop.
			}
		} // strand signal check.

		delete(cur_entry_length_str);
		delete [] cur_read;
	} // file reading loop.
	fclose(f_sorted_reads);

	// Get the length of the data.
	l_data = l_buffer;
	while(l_data > 0)
	{
		if(signal_profile_buffer != NULL && signal_profile_buffer[l_data] > 0.0)
		{
			break;
		}
		else if(forward_strand_signal != NULL && forward_strand_signal[l_data] > 0.0)
		{
			break;
		}
		else
		{
			l_data--;
		}
	} // l_data setting loop.
}

int get_l_signal_per_reads(char* reads_fp, int l_ext_tag)
{
	char strand_char = 0;
	//char cur_fragment[10000];
	char mapping_map_str[10000];
	int chr_index;
	int n_processed_reads = 0;
	int l_profile = 0;
	FILE* f_sorted_reads = open_f(reads_fp, "r");
	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
		}	

		n_processed_reads++;

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int cur_l_profile = chr_index + 1000 + l_ext_tag;
		if(cur_l_profile > l_profile)
		{
			l_profile = cur_l_profile;
		}
	} // file loading loop.

	fclose(f_sorted_reads);

	return(l_profile);
} // get_l_signal_per_reads option.

// In order to ignore extended tag length, set it to non-positive value.
void buffer_per_nucleotide_profile_no_buffer(const char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data)
{
	for(int i = 0; i <= l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initilize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	// Non-positive tag extension lengths are ignored and falls back to using the length in the CIGAR string entry.
	if(l_extended_tag <= 0)
	{
		fprintf(stderr, "Ignoring tag extension.\n");
	}

	char strand_char = 0;
	//char cur_fragment[100000];
	char mapping_map_str[100000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = open_f(sorted_read_fp, "r");

	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
}
		}	

		n_processed_reads++;

		// We need a check on the current read line to make sure what we have is a valid: index must be strictly numbers; mapping map string is validated
		// below.
		char chr_index_str[1000];

		if(sscanf(cur_read, "%s %c %s", mapping_map_str, &strand_char, chr_index_str) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		// Validate the string read as chromosome index.
		int char_i = 0;
		while(chr_index_str[char_i] != 0)
		{
			bool char_is_a_num = (chr_index_str[char_i] >= '0' && chr_index_str[char_i] <= '9');
			if(!char_is_a_num)
			{
				fprintf(stderr, "Chromosome index must be a number: %s\n", cur_read);
				exit(0);
			}

			char_i++;
		} // char_i loop.

		// This may be pointing to an error in read preprocessing.
		chr_index = atoi(chr_index_str);
		if(chr_index > l_buffer)
		{
			fprintf(stderr, "%s: The read mapped out of buffer.\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		int right_most_match_pos = 0;
		int left_most_match_pos = 1000*1000*1000;

		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			if(is_matching)
			{
				// Should we ignore the tag extension?
				if(l_extended_tag <= 0)
				{
					// Increase the height: The indexing is already 1 based, there is no conversion needed.
					if(signal_profile_buffer != NULL)
					{
						for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
						{
							signal_profile_buffer[i_cur]++;
						} // i_cur loop.
					}

					// Update the strand signals, if requested.
					if(forward_strand_signal != NULL)
					{
						if(strand_char == 'F')
						{
							// Update the forward strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								forward_strand_signal[i_cur]++;
							} // i_cur loop.
						}
						else
						{
							// Update the reverse strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								reverse_strand_signal[i_cur]++;
							} // i_cur loop.
						}
					} // strand signal check.
				}
				else // tag extension validation check.
				{
					left_most_match_pos = MIN(left_most_match_pos, chr_index);
					right_most_match_pos = (chr_index + l_cur_entry - 1);
				} // extension length check.
			} // match block check.

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// If tag extension is requested, utilize the leftmost and rightmost matching position for the mapped read.
		if(l_extended_tag > 0)
		{
			int ext_start = 0;
			int ext_end = 0;
			if(strand_char == 'F')
			{
				ext_start = left_most_match_pos;
			}
			else
			{
				//ext_start = (chr_index + l_cur_entry - 1) - (l_extended_tag - 1);
				ext_start = right_most_match_pos - (l_extended_tag - 1);
			}

			// Check for negative starts.
			if(ext_start < 0)
			{
				ext_start = 1;
			}

			ext_end = ext_start + l_extended_tag - 1;

			// Update profiles for the strands.
			if(forward_strand_signal != NULL)
			{
				if(strand_char == 'F')
				{
					// Update the forward strand signal.
					for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
					{
						forward_strand_signal[i_cur]++;
					} // i_cur loop.
				}
				else
				{
					// Update the reverse strand signal.
					for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
					{
						reverse_strand_signal[i_cur]++;
					} // i_cur loop.
				}
			} // strand signal check.

			// Update the profile.
			if(signal_profile_buffer != NULL)
			{
				for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
				{
					signal_profile_buffer[i_cur]++;
				} // i_cur loop.
			} // profile existence check.
		} // signal update check.

		delete [] cur_read;
	} // file reading loop.

	close_f(f_sorted_reads, sorted_read_fp);

	// Get the length of the data.
	l_data = l_buffer;
	while(l_data > 0)
	{
		if(signal_profile_buffer != NULL && signal_profile_buffer[l_data] > 0.0)
		{
			break;
		}
		else if(forward_strand_signal != NULL && forward_strand_signal[l_data] > 0.0)
		{
			break;
		}
		else
		{
			l_data--;
		}
	} // l_data setting loop.
}

// In order to ignore extended tag length, set it to non-positive value.
void buffer_per_nucleotide_profile_no_buffer_ORIGINAL_NO_SPLICE_PROCESSING(char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data)
{
	for(int i = 0; i < l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initilize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	// Non-positive tag extension lengths are ignored and falls back to using the length in the CIGAR string entry.
	if(l_extended_tag <= 0)
	{
		fprintf(stderr, "Ignoring tag extension.\n");
	}

	char strand_char = 0;
	//char cur_fragment[100000];
	char mapping_map_str[100000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = NULL;
	if(t_string::compare_strings(sorted_read_fp, "stdin"))
	{
		f_sorted_reads = stdin;
	}
	else
	{
		f_sorted_reads = open_f(sorted_read_fp, "r");
	}

	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
		}	

		n_processed_reads++;

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
		bool splicing_tag_extension_check_pass = true;
		if(is_read_spliced && l_extended_tag > 0)
		{
			fprintf(stderr, "Spliced read with tag extension: %s\r", mapping_map_str);
			splicing_tag_extension_check_pass = false;
		}

		while(splicing_tag_extension_check_pass &&
			mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			if(is_matching)
			{
				// Should we ignore the tag extension?
				if(l_extended_tag <= 0)
				{
					// Increase the height: The indexing is already 1 based, there is no conversion needed.
					if(signal_profile_buffer != NULL)
					{
						for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
						{
							signal_profile_buffer[i_cur]++;
						} // i_cur loop.
					}

					// Update the strand signals, if requested.
					if(forward_strand_signal != NULL)
					{
						if(strand_char == 'F')
						{
							// Update the forward strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								forward_strand_signal[i_cur]++;
							} // i_cur loop.
						}
						else
						{
							// Update the reverse strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								reverse_strand_signal[i_cur]++;
							} // i_cur loop.
						}
					} // strand signal check.
				}
				else // tag extension validation check.
				{
					// Extend the read, update the profile.
					int ext_start = 0;
					int ext_end = 0;
					if(strand_char == 'F')
					{
						ext_start = chr_index;
					}
					else
					{
						ext_start = (chr_index + l_cur_entry - 1) - (l_extended_tag - 1);
					}

					// Check for negative starts.
					if(ext_start < 0)
					{
						ext_start = 1;
					}

					ext_end = ext_start + l_extended_tag - 1;

					// Update profiles for the strands.
					if(forward_strand_signal != NULL)
					{
						if(strand_char == 'F')
						{
							// Update the forward strand signal.
							for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
							{
								forward_strand_signal[i_cur]++;
							} // i_cur loop.
						}
						else
						{
							// Update the reverse strand signal.
							for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
							{
								reverse_strand_signal[i_cur]++;
							} // i_cur loop.
						}
					} // strand signal check.

					// Update the profile.
					if(signal_profile_buffer != NULL)
					{
						for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
						{
							signal_profile_buffer[i_cur]++;
						} // i_cur loop.
					}
				} // extension length check.
			}

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		delete(cur_entry_length_str);
		delete [] cur_read;
	} // file reading loop.

	if(t_string::compare_strings(sorted_read_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_sorted_reads);
	}

	// Get the length of the data.
	l_data = l_buffer;
	while(l_data > 0)
	{
		if(signal_profile_buffer != NULL && signal_profile_buffer[l_data] > 0.0)
		{
			break;
		}
		else if(forward_strand_signal != NULL && forward_strand_signal[l_data] > 0.0)
		{
			break;
		}
		else
		{
			l_data--;
		}
	} // l_data setting loop.
}

void delete_fragments(vector<t_mapped_fragment*>* fragment_list)
{
	for(int i_frag = 0; i_frag < (int)fragment_list->size(); i_frag++)
	{
		delete(fragment_list->at(i_frag));
	}

	fragment_list->clear();
	delete(fragment_list);
}

void delete_fragments(t_mapped_fragment** fragment_list)
{
	int i_frag = 0;
        while(fragment_list[i_frag] != NULL)
        {
                free(fragment_list[i_frag]);
		i_frag++;
        }

        delete[] fragment_list;
}

void load_fragments(char* mapped_reads_fp, 
	vector<t_mapped_fragment*>* fore_strand_frags, vector<t_mapped_fragment*>* rev_strand_frags, 
	int max_n_pcr_amplified_reads)
{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
	fprintf(stderr, "Loading fragments from %s.\n", mapped_reads_fp);
}

	vector<t_mapped_read*>* fore_reads = new vector<t_mapped_read*>();
	vector<t_mapped_read*>* rev_reads = new vector<t_mapped_read*>();

	// Add the fragments: If the fragmens do not exist, the algorithm returns.
	//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
	load_reads(mapped_reads_fp, fore_reads, rev_reads, max_n_pcr_amplified_reads);
	get_mapped_fragments_per_mapped_reads(fore_reads, fore_strand_frags);
	get_mapped_fragments_per_mapped_reads(rev_reads, rev_strand_frags);

	// Delete the memory for the mapped reads, that is not needed any more.
	for(int i_f_r = 0; i_f_r < (int)fore_reads->size(); i_f_r++)
	{
		delete_mapped_read(fore_reads->at(i_f_r));
	} // i_f_r loop.
	delete fore_reads;

	for(int i_r_r = 0; i_r_r < (int)rev_reads->size(); i_r_r++)
	{
		delete_mapped_read(rev_reads->at(i_r_r));
	} // i_f_r loop.
	delete rev_reads;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
	fprintf(stderr, "Loaded %ld forward and %ld reverse fragments.\n", fore_strand_frags->size(), rev_strand_frags->size());
}
}

bool check_fragment_quality(char* fragment, char* chr_subseq, char* quality_str)
{
        int mismatches = 0;
        for(int i = 0; i < (int)strlen(fragment); i++)
        {
                char upper_subseq = toupper(chr_subseq[i]);
                char upper_frag = toupper(fragment[i]);

                if((upper_subseq != upper_frag) &&
                        (upper_subseq != 'N') &&
                        (upper_frag != 'N'))
                {
                        mismatches++;
                        //printf("Fragments are not matching for %s at %d (%s). Enriched fragment is:\n%s\n%s\n", chr_fps->at(i_chr), chr_index, quality_str, cur_seq, chr_subseq);
                        //getc(stdin);
                        //exit(0);
                }
        } // i loop.

        if(mismatches == 0)
        {
		return true;
/*
                if(strcmp(quality_str, "U0") == 0)
                {
                        return(true);
                }
                else
                {
                        return(false);
                }
*/
        }
        else if(mismatches == 1)
        {
                if(strcmp(quality_str, "U1") == 0)
                {
                        return(true);
                }
                else
                {
                        return(false);
                }
        }
        else if(mismatches == 2)
        {
                if(strcmp(quality_str, "U2") == 0)
                {
                        return(true);
                }
                else
                {
                        return(false);
                }
        }
        else
        {
		return(false);
        }
}

bool sort_mapped_fragments_per_3p(t_mapped_fragment* frag1, t_mapped_fragment* frag2)
{
	int frag1_3p = fragment_3p_accessor(frag1);
	int frag2_3p = fragment_3p_accessor(frag2);

	return(frag1_3p < frag2_3p);
}

bool sort_mapped_fragments(t_mapped_fragment* frag1, t_mapped_fragment* frag2)
{
	if(frag1->base_index < frag2->base_index)
	{
		return(frag1->base_index < frag2->base_index);
	}
	else if(frag1->base_index == frag2->base_index)
	{
		if((frag1->base_index + frag1->sequenced_fragment_length) < (frag2->base_index + frag2->sequenced_fragment_length))
		{
			return(true);
		}
		else
		{
			return(false);
		}
	}
	else
	{
		return(false);
	}
}

int fragment_5p_accessor(void* obj_ptr)
{
	t_mapped_fragment* frag_obj_ptr = (t_mapped_fragment*)obj_ptr;

	if(frag_obj_ptr->strand_char == 'F')
	{
		return(frag_obj_ptr->base_index);		
	}
	else if(frag_obj_ptr->strand_char == 'R')
	{
		return(frag_obj_ptr->base_index);
	}
	else
	{
		fprintf(stderr, "The strand char for fragment object is %c @ %s(%d)\n", frag_obj_ptr->strand_char, __FILE__, __LINE__);
		exit(0);
		return(-1);
	}
}

int fragment_3p_accessor(void* obj_ptr)
{
	t_mapped_fragment* frag_obj_ptr = (t_mapped_fragment*)obj_ptr;

	if(frag_obj_ptr->strand_char == 'F')
	{
		return(frag_obj_ptr->base_index+frag_obj_ptr->sequenced_fragment_length-1);
	}
	else if(frag_obj_ptr->strand_char == 'R')
	{
		return(frag_obj_ptr->base_index+frag_obj_ptr->sequenced_fragment_length-1);
	}
	else
	{
		fprintf(stderr, "The strand char for fragment object is %c @ %s(%d)\n", frag_obj_ptr->strand_char, __FILE__, __LINE__);
		exit(0);
		return(-1);
	}
}

// Count the statistics of mapped fragments for a given regions using binary search.
void get_read_statistics_per_region(vector<t_mapped_fragment*>* fragments, 
									vector<int>* all_fragment_3p_posns, 
									t_annot_region* region, 
									double& n_mapped_reads, 
									double& n_mapped_nucs)
{
	n_mapped_nucs = 0;
	n_mapped_reads = 0;

	// For the given region, count the reads with a binary search. Make sure that the fragments are sorted.
	int i_leftmost_read = locate_posn_per_sorted_posn_list(region->start, all_fragment_3p_posns, 0, all_fragment_3p_posns->size()-1);
	//int i_leftmost_read = locate_posn_per_sorted_obj_list(region->start, (vector<void*>*)all_fragment_3p_posns, 0, all_fragment_3p_posns->size()-1, fragment_3p_accessor);

	//fprintf(stderr, "Found leftmost read (%d-%d) @ %d. read for region %s:%d-%d\n", fragments->at(i_leftmost_read)->base_index, fragments->at(i_leftmost_read)->base_index + fragments->at(i_leftmost_read)->sequenced_fragment_length, i_leftmost_read, region->chrom, region->start, region->end);

	// Starting from the leftmost read, count the overlapping reads.
	while(i_leftmost_read < (int)fragments->size() && 
		fragments->at(i_leftmost_read)->base_index <= region->end)
	{
		// Get the overlap of the current read with the current region.
		int ol_start = MAX(fragments->at(i_leftmost_read)->base_index, region->start);
		int ol_end = MIN(fragments->at(i_leftmost_read)->base_index + fragments->at(i_leftmost_read)->sequenced_fragment_length-1, region->end);

		if(ol_end >= ol_start)
		{
			n_mapped_reads++;
			n_mapped_nucs += (ol_end - ol_start + 1);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "%s:%d-%d: Read @ %d, Overlap: %d-%d Totals: %d, %d\n", region->chrom, region->start, region->end, fragments->at(i_leftmost_read)->base_index, ol_start, ol_end, (int)n_mapped_nucs, (int)n_mapped_reads);
		}

		i_leftmost_read++;
	} // i_leftmost_read loop.
}

void load_reads_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
	vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
	int max_n_pcr_amplified_reads)
{
	char cur_dir_chr_ids_fp[1000];
	sprintf(cur_dir_chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	vector<char*>* cur_dir_chr_ids = buffer_file(cur_dir_chr_ids_fp);

	if(cur_dir_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", mapped_reads_dir);
		getc(stdin);
		return;
	}

	// Load the fragment in this chromosome from all the directories.
	for(int i_chr = 0; i_chr < (int)cur_dir_chr_ids->size(); i_chr++)
	{
		//fore_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
		//rev_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());

		chr_ids->push_back(t_string::copy_me_str(cur_dir_chr_ids->at(i_chr)));

		char cur_dir_cur_chr_reads_fp[1000];
		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, cur_dir_chr_ids->at(i_chr));

		vector<t_mapped_read*>* cur_chr_fore_reads = new vector<t_mapped_read*>();
		vector<t_mapped_read*>* cur_chr_rev_reads = new vector<t_mapped_read*>();

		// Add the fragments: If the fragmens do not exist, the algorithm returns.
		//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
		load_reads(cur_dir_cur_chr_reads_fp, cur_chr_fore_reads, cur_chr_rev_reads, max_n_pcr_amplified_reads);

		// Add the loaded fragments.
		fore_strand_reads_per_chr->push_back(cur_chr_fore_reads);
		rev_strand_reads_per_chr->push_back(cur_chr_rev_reads);
	} // i_chr loop.
}

void load_fragments_per_dir(char* mapped_reads_dir, 
	vector<char*>* chr_ids, 
	vector<vector<t_mapped_fragment*>*>* fore_strand_fragments_per_chr, 
	vector<vector<t_mapped_fragment*>*>* rev_strand_fragments_per_chr, 
	int max_n_pcr_amplified_reads)
{
	char cur_dir_chr_ids_fp[1000];
	sprintf(cur_dir_chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
	vector<char*>* cur_dir_chr_ids = buffer_file(cur_dir_chr_ids_fp);

	if(cur_dir_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", mapped_reads_dir);
		getc(stdin);
		return;
	}

	// Load the fragment in this chromosome from all the directories.
	for(int i_chr = 0; i_chr < (int)cur_dir_chr_ids->size(); i_chr++)
	{
		//fore_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
		//rev_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());

		chr_ids->push_back(t_string::copy_me_str(cur_dir_chr_ids->at(i_chr)));

		char cur_dir_cur_chr_reads_fp[1000];
		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, cur_dir_chr_ids->at(i_chr));

		vector<t_mapped_read*>* cur_chr_fore_reads = new vector<t_mapped_read*>();
		vector<t_mapped_read*>* cur_chr_rev_reads = new vector<t_mapped_read*>();

		// Add the fragments: If the fragmens do not exist, the algorithm returns.
		//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
		load_reads(cur_dir_cur_chr_reads_fp, cur_chr_fore_reads, cur_chr_rev_reads, max_n_pcr_amplified_reads);

		vector<t_mapped_fragment*>* cur_chr_fore_frags = new vector<t_mapped_fragment*>(); 
		get_mapped_fragments_per_mapped_reads(cur_chr_fore_reads, cur_chr_fore_frags);

		vector<t_mapped_fragment*>* cur_chr_rev_frags = new vector<t_mapped_fragment*>(); 
		get_mapped_fragments_per_mapped_reads(cur_chr_rev_reads, cur_chr_rev_frags);

		// Add the loaded fragments.
		fore_strand_fragments_per_chr->push_back(cur_chr_fore_frags);
		rev_strand_fragments_per_chr->push_back(cur_chr_rev_frags);

		// Delete the memory for the mapped reads, that is not needed any more.
		delete_mapped_reads(cur_chr_fore_reads);
		delete_mapped_reads(cur_chr_rev_reads);
	} // i_chr loop.
}

vector<t_mapped_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_mapped_fragment*>* fore_frag_list, vector<t_mapped_fragment*>* rev_frag_list, int enrichment_mapped_fragment_length)
{
	vector<t_mapped_fragment*>* combined_frags = new vector<t_mapped_fragment*>();

	if(fore_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < (int)fore_frag_list->size(); i_frag++)
		{
			t_mapped_fragment* new_fragment_node = (t_mapped_fragment*)malloc(sizeof(t_mapped_fragment));
			new_fragment_node->base_index = fore_frag_list->at(i_frag)->base_index;

			// If the enrichment length is greater than the actual sequence tag length, update the sequenced length.
			if(enrichment_mapped_fragment_length > fore_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				new_fragment_node->sequenced_fragment_length = enrichment_mapped_fragment_length;
			}
			else
			{
				new_fragment_node->sequenced_fragment_length = fore_frag_list->at(i_frag)->sequenced_fragment_length;
			}


			new_fragment_node->strand_char = fore_frag_list->at(i_frag)->strand_char;

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// For the reverse fragment list, we need to also set the base index since it may be move further to left.
	if(rev_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < (int)rev_frag_list->size(); i_frag++)
		{
			// If this enrichment fragment length is set to 0, this means that there is no fragment extension.
			t_mapped_fragment* new_fragment_node = (t_mapped_fragment*)malloc(sizeof(t_mapped_fragment));

			// If the enrichment length is greater than the actual sequence tag length, update the sequenced length.
			if(enrichment_mapped_fragment_length > rev_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				// The enrichment fragment length is larger than the sequenced length, translate the fragment base.
				new_fragment_node->base_index = (rev_frag_list->at(i_frag)->base_index + rev_frag_list->at(i_frag)->sequenced_fragment_length - enrichment_mapped_fragment_length);
				new_fragment_node->sequenced_fragment_length = enrichment_mapped_fragment_length;
			}
			else
			{
				// The enrichment fragment length is smaller than the sequenced fragment length, do not do any translation of the fragment base.
				new_fragment_node->base_index = rev_frag_list->at(i_frag)->base_index;
				new_fragment_node->sequenced_fragment_length = rev_frag_list->at(i_frag)->sequenced_fragment_length;
			}

			new_fragment_node->strand_char = 'F'; // Change the place of these to forward strand.

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// Sort one last time.
	sort(combined_frags->begin(), combined_frags->end(), sort_mapped_fragments);

	return(combined_frags);
}

int* get_n_reads_per_window(int n_wins, vector<t_mapped_fragment*>* frag_list)
{
	int* n_reads_per_window = new int[n_wins + 2];
	for(int i_win = 0; i_win < n_wins; i_win++)
	{
		n_reads_per_window[i_win] = 0;
	} // i_win loop

	for(int i_frag = 0; i_frag < (int)frag_list->size(); i_frag++)
	{
		int cur_i_win = frag_list->at(i_frag)->base_index / MEG_BASE;

		n_reads_per_window[cur_i_win]++;
	} // i_nuc loop.

	return(n_reads_per_window);
}

enum{VAL, TYPE};
bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced)
{
	int i = 0;

	is_read_spliced = false;

	int state = VAL;
	while(mapping_map_str[i] != 0)
	{		
		if(state == VAL)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 0);
			// MIDNSHPX=
			if(mapping_map_str[i] == 'M' ||
				mapping_map_str[i] == 'I' ||
				mapping_map_str[i] == 'D' ||
				mapping_map_str[i] == 'N' ||
				mapping_map_str[i] == 'S' ||		
				mapping_map_str[i] == 'H' ||
				mapping_map_str[i] == 'P' ||
				mapping_map_str[i] == 'X' ||
				mapping_map_str[i] == '=')
			{
				state = TYPE;

				//if(mapping_map_str[i] != 'M')
				// If the state is N, we assume that this is a spliced read: Page 7, CIGAR definition @ https://samtools.github.io/hts-specs/SAMv1.pdf 
				if (mapping_map_str[i] == 'N')
				{
					is_read_spliced = true;
				}
			}
			else if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				// State is still VAL.
			}
			else
			{
				return(false);
			}
		}
		else if(state == TYPE)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 1);
			// A number is expected.
			if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				state = VAL;
			}
			else
			{
				return(false);
			}
		}

		// Move to next character.
		i++;
	}

	return(true);
}

void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										//t_string* cur_entry_length_str,
										int& l_cur_entry,
										char& entry_type_char)
{	
	// Clean the length string.
	//cur_entry_length_str->empty();
	l_cur_entry = 0;

	// Get the next entry in the cigar string.
	while(mapping_map_str[i_mapp_map] != 0)
	{
		if(mapping_map_str[i_mapp_map] < '0' || mapping_map_str[i_mapp_map] > '9')
		{
			break;
		}
		//cur_entry_length_str->concat_char(mapping_map_str[i_mapp_map]);
		l_cur_entry = l_cur_entry*10 + (int)(mapping_map_str[i_mapp_map]-'0');
		i_mapp_map++;
	}

	is_matching = false;
	if(mapping_map_str[i_mapp_map] == 'M')
	{
		//fprintf(stderr, "Adding matching length of %d\n", l_cur_entry);
		is_matching = true;
	}
	else
	{
		//fprintf(stderr, "Adding some other length of %d\n", l_cur_entry);
	}	

	entry_type_char = mapping_map_str[i_mapp_map];

	// Move over the current entry identifier.
	i_mapp_map++;
}

void load_PE_reads(char* sorted_first_mapped_reads_fp, 
	char* sorted_last_mapped_reads_fp, 
	vector<t_mapped_PE_read*>* first_mate_reads, 
	vector<t_mapped_PE_read*>* last_mate_reads, 
	int max_n_pcr_amplified_reads)
{
	if(!check_file(sorted_first_mapped_reads_fp))
	{
		printf("Could not open mapped reads file %s @ %s(%d).\n", sorted_first_mapped_reads_fp, __FILE__, __LINE__);
		return;
	}

	if(!check_file(sorted_last_mapped_reads_fp))
	{
		printf("Could not open mapped reads file %s @ %s(%d).\n", sorted_first_mapped_reads_fp, __FILE__, __LINE__);
		return;
	}

	// Start loading the PE reads: For each read on the first mate file, find the corresponding id'd read.
}

/*
Load reads from a preprocesed read file.
*/
void load_reads(char* mapped_reads_fp, 
	vector<t_mapped_read*>* pruned_fore_strand_reads, vector<t_mapped_read*>* pruned_rev_strand_reads, 
	int max_n_pcr_amplified_reads)
{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
	printf("Loading mapped-reads from %s.\n", mapped_reads_fp);

	if(!check_file(mapped_reads_fp))
	{
		printf("Could not open mapped reads file %s @ %s(%d).\n", mapped_reads_fp, __FILE__, __LINE__);
		return;
	}

	vector<t_mapped_read*>* fore_strand_reads = new vector<t_mapped_read*>();
	vector<t_mapped_read*>* rev_strand_reads = new vector<t_mapped_read*>(); 

	//char cur_fragment[10000];
	char mapping_map_str[10000];
	char strand_char;
	int chr_index;

	// Buffer the whole file.
	FILE* f_mrf = open_f(mapped_reads_fp, "r");

	// Read and validate the mapped reads in the file.
	//while(fscanf(f_mapped_reads, "%s %s %c %d", cur_fragment, quality_str, &strand_char, &chr_index) == 4)
	while(1)
	{
		char* cur_line = getline(f_mrf);

		if(cur_line == NULL)
		{
			break;
		}

		if(sscanf(cur_line, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse fragment line: %s\n", cur_line);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		//fprintf(stderr, "Processing cigar string: %s\n", quality_str);
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		// Allocate and initialize the new mapped read.
		t_mapped_read* new_mapped_read = new t_mapped_read();
		new_mapped_read->mapping_str = t_string::copy_me_str(mapping_map_str);
		new_mapped_read->strand = strand_char;
		new_mapped_read->span = 0;
		int left_posn = chr_index;

		// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
				new_mapped_read->span += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Set the base_index.
		new_mapped_read->base_index = left_posn;	

		new_mapped_read->strand = new_mapped_read->strand;

		// Add to the lists.
		if(new_mapped_read->strand == 'F')
		{
			fore_strand_reads->push_back(new_mapped_read);
		}
		else if(new_mapped_read->strand == 'R')
		{
			rev_strand_reads->push_back(new_mapped_read);
		}

		delete(cur_entry_length_str);
		delete [] cur_line;
	} // current fragment data reading loop.

	// Free file buffer.
	close_f(f_mrf, mapped_reads_fp);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
	printf("Loaded %ld fragments on forward strand.\n", fore_strand_reads->size());
    printf("Loaded %ld fragments on reverse strand.\n", rev_strand_reads->size());
}

	// Prune the reads.
	vector<t_mapped_read*>* all_reads = new vector<t_mapped_read*>();
	all_reads->insert(all_reads->end(), fore_strand_reads->begin(), fore_strand_reads->end());
	all_reads->insert(all_reads->end(), rev_strand_reads->begin(), rev_strand_reads->end());
	
	// Following gets rid of the pruned reads and returns only the pruned reads without re-allocating the reads/split-mapped-fragments.
	prune_reads(all_reads, max_n_pcr_amplified_reads, 
				pruned_fore_strand_reads, 
				pruned_rev_strand_reads);

	// Can clean the memory for all reads. Note that the pruned reads are cleaned up above, unpruned ones are put in the reads per strand lists.
	delete(fore_strand_reads);
	delete(rev_strand_reads);
	delete(all_reads);
}

// This is a generic iterator over the processed reads file, useful for doing online processing of the reads files, for example, when they are too large to load into memory.
void preprocessed_read_file_iterator(char* mapped_reads_fp,
	void (per_read_callback)(char*, char, int, void*), 
	void (per_fragment_callback)(char*, char, int, void*),
	void* per_read_callback_param,
	void* per_fragment_callback_param)
{
	printf("Loading mapped-reads from %s.\n", mapped_reads_fp);

	FILE* f_mrf = open_f(mapped_reads_fp, "r");

	//char cur_fragment[10000];
	char mapping_map_str[10000];
	char strand_char;
	int chr_index;

	// Read and validate the mapped reads in the file.
	//while(fscanf(f_mapped_reads, "%s %s %c %d", cur_fragment, quality_str, &strand_char, &chr_index) == 4)
	while(1)
	{
		char* cur_line = getline(f_mrf);

		if(cur_line == NULL)
		{
			break;
		}

		if(sscanf(cur_line, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse fragment line: %s\n", cur_line);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		//fprintf(stderr, "Processing cigar string: %s\n", quality_str);
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		//int left_posn = chr_index;

		// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
			if(is_matching && per_fragment_callback != NULL)
			{
				// Call the fragment callback.
				per_fragment_callback(mapping_map_str, strand_char, chr_index, per_fragment_callback_param);
			}

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Call the read callback.
		per_read_callback(mapping_map_str, strand_char, chr_index, per_read_callback_param);

		//fprintf(stderr, "%s: %s, %d (%d), %c\n", cur_line, new_mapped_read->mapping_str, new_mapped_read->base_index, new_mapped_read->span, new_mapped_read->strand);
		//getc(stdin);

		delete(cur_entry_length_str);
		delete [] cur_line;
	} // curent fragment data reading loop.

	close_f(f_mrf, mapped_reads_fp);
}

void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments)
{
	int i_mapp_map = 0;
	t_string* cur_entry_length_str = new t_string();
	bool is_matching = false;
	char entry_type_char;
	char strand_char = mapped_read->strand;
	//int chr_index = (mapped_read->strand=='F')?(mapped_read->base_index):(mapped_read->base_index-mapped_read->span+1);
	int chr_index = (mapped_read->base_index);
	char* mapping_map_str = mapped_read->mapping_str;

	//fprintf(stderr, "Processing cigar string: %s (%d, %c)\n", mapping_map_str, chr_index, strand_char);
	bool is_read_spliced = false;
	bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

	// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
	while(mapping_map_str_valid && 
		mapping_map_str[i_mapp_map] != 0)
	{
		int l_cur_entry = 0;
		get_next_entry_per_mapp_map_string(mapping_map_str,
											i_mapp_map, 
											is_matching,
											l_cur_entry,
											entry_type_char);

		// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.		

		//int l_fragment = strlen(cur_fragment);
		if(is_matching)
		{
			//int l_fragment = get_l_fragment_per_cigar(quality_str);
			if(strand_char == 'F')
			{
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;
		
				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
			}
			else if(strand_char == 'R')
			{
				// Allocate and initialize a fragment and add it to the reverse strand fragment list.			
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				//new_fragment->base_index = chr_index + l_cur_entry - 1;
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;

				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
				//rev_strand_frags->push_back(new_fragment);
			} // reverse strand check.
		} // maching check.

		// Update the base for the current entry.
		// Must check whether to update the mapping posn: Update only for D and M entries.
		//if(entry_type_char == 'D' || 
		//	entry_type_char == 'M' ||
		//	entry_type_char == 'N' ||
		//	entry_type_char == 'H')
		if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
		{
			chr_index += l_cur_entry;
		}
	} // mapping map string processing loop.

	delete cur_entry_length_str;
}

// Following should be current with SAM specs.
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char)
{
	if(entry_char == '=' ||
		entry_char == 'X' ||
		entry_char == 'S' ||
		entry_char == 'M' ||
		entry_char == 'I')
	{
		return(true);
	}

	return(false);
}

// Following should be current with SAM specs.
bool check_genome_index_update_per_CIGAR_entry(char entry_char)
{
	if(entry_char == 'D' || 
		entry_char == 'M' ||
		entry_char == 'N' ||
		entry_char == '=' ||
		entry_char == 'X')
	{
		return(true);
	}

	return(false);
}

void exclude_reads_per_regions_per_chr(char* read_chr_id, 
	vector<t_mapped_read*>* cur_chr_reads,
	vector<t_mapped_read*>* no_overlap_reads,
	vector<t_annot_region*>* regions_2_exclude)
{
	if(regions_2_exclude == NULL)
	{
		for(int i_r = 0; i_r < (int)cur_chr_reads->size(); i_r++)
		{
			no_overlap_reads->push_back(cur_chr_reads->at(i_r));
		} // i_r loop.

		return;
	}

	// Restructure the regions.
	t_restr_annot_region_list* restructured_regions = restructure_annot_regions(regions_2_exclude);
	
	int n_total_excluded_reads = 0;
	
	int n_excluded_reads = 0;

	int i_reg_chr = t_string::get_i_str(restructured_regions->chr_ids, read_chr_id);

	if(i_reg_chr < (int)restructured_regions->chr_ids->size())
	{
		fprintf(stderr, "Excluding the reads in %s\n", restructured_regions->chr_ids->at(i_reg_chr));			

		for(int i_read = 0; i_read < (int)cur_chr_reads->size(); i_read++)
		{
			int cur_left_base_posn = (cur_chr_reads->at(i_read)->base_index);
			int cur_right_base_posn = (cur_chr_reads->at(i_read)->base_index + cur_chr_reads->at(i_read)->span - 1);

			int left_reg_i = locate_posn_per_sorted_obj_list(cur_left_base_posn, (vector<void*>*)(restructured_regions->regions_per_chrom[i_reg_chr]), 0, restructured_regions->regions_per_chrom[i_reg_chr]->size()-1, region_5p_accessor);

			bool cur_read_overlaps = false;
			while(left_reg_i > 0 &&
				restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end > cur_left_base_posn)
			{
				left_reg_i--;
			}

			while(left_reg_i < (int)restructured_regions->regions_per_chrom[i_reg_chr]->size() &&
				restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start < cur_right_base_posn)
			{
				int ol_start = MAX(restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start, cur_left_base_posn);
				int ol_end = MIN(restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end, cur_right_base_posn);

				if(ol_end >= ol_start)
				{
					cur_read_overlaps = true;
					break;
				}

				left_reg_i++;
			} // go over the region for the read.

			if(cur_read_overlaps)
			{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				if(cur_chr_reads->at(i_read)->strand == 'R')
				{
					fprintf(stderr, "Read: %s:%d(%d)(%c) overlaps with %s:%d-%d\n", read_chr_id, cur_chr_reads->at(i_read)->base_index, cur_chr_reads->at(i_read)->span, cur_chr_reads->at(i_read)->strand, 
						restructured_regions->chr_ids->at(i_reg_chr), restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start, restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end);
					getc(stdin);
				}
}

				// Update # of excluded reads.
				n_excluded_reads++;
			}
			else
			{
				no_overlap_reads->push_back(cur_chr_reads->at(i_read));
			}
		} // i_read loop.
	} // i_reg_chr search check.
	else
	{
		fprintf(stderr, "Skipping the regions in %s\n", read_chr_id);
		for(int i_read = 0; i_read < (int)cur_chr_reads->size(); i_read++)
		{
			no_overlap_reads->push_back(cur_chr_reads->at(i_read));
		} // i_read loop.
	} // read pile chromosome check.

	fprintf(stderr, "Excluded %d reads.\n", n_excluded_reads);

	n_total_excluded_reads += n_excluded_reads;

	delete_restructured_annot_regions(restructured_regions);

	fprintf(stderr, "Excluded %d reads in total.\n", n_total_excluded_reads);
}

void exclude_reads_per_regions(vector<char*>* read_chr_ids, 
	vector<vector<t_mapped_read*>*>* reads_per_chrs, 
	vector<vector<t_mapped_read*>*>* no_overlap_reads_per_chrs,
	vector<t_annot_region*>* regions_2_exclude)
{
	for(int i_read_chr = 0; i_read_chr < (int)read_chr_ids->size(); i_read_chr++)
	{
		vector<t_mapped_read*>* cur_chr_no_overlap_reads = new vector<t_mapped_read*>();
		no_overlap_reads_per_chrs->push_back(cur_chr_no_overlap_reads);

		exclude_reads_per_regions_per_chr(read_chr_ids->at(i_read_chr), 
			reads_per_chrs->at(i_read_chr),
			cur_chr_no_overlap_reads,
			regions_2_exclude);
	} // i_read_chr loop.
}

// Merge all the fragments of all the reads.
void get_mapped_fragments_per_mapped_reads(vector<t_mapped_read*>* mapped_reads, vector<t_mapped_fragment*>* mapped_fragments)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		add_mapped_fragments_per_mapped_read(mapped_reads->at(i_r), mapped_fragments);
	} // i_r loop.
}

void delete_mapped_reads(vector<t_mapped_read*>* mapped_reads)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		delete_mapped_read(mapped_reads->at(i_r));
	} // i_r loop.

	delete(mapped_reads);
}

void delete_mapped_read(t_mapped_read* mapped_read)
{
	delete [] mapped_read->mapping_str;
	delete(mapped_read);
}

/*
Prune reads:
Note that the trick here is to deal with the reads that have exact same pattern of mapping. We do not care about the
strand, this should be taken care of before the function is called.

Note that the pruning must be done at the read level, not at the fragment level.
*/
void prune_reads(vector<t_mapped_read*>* mapped_reads, int n_max_reps_per_posn, 
	vector<t_mapped_read*>* pruned_forward_reads, 
	vector<t_mapped_read*>* pruned_reverse_reads)
{
	// If the pruning is not requested, return all the reads.
	if(n_max_reps_per_posn == 0)
	{
		fprintf(stderr, "Skipping pruning.\n");
		for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		} // i_r loop.

		return;
	}

	// Sort the mapped reads with respect to their 5' posn.
	sort(mapped_reads->begin(), mapped_reads->end(), sort_mapped_reads_per_5p);

    // First get rid of the extra fragments on forward strand.
	int* rep_cnts = new int[mapped_reads->size() + 2];
	memset(rep_cnts, 0, sizeof(int) * (mapped_reads->size() + 1));

	int prev_read_left_posn = 0;
    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		int cur_read_left_posn = (mapped_reads->at(i_r)->base_index);
		
		if(i_r > 0 &&
			prev_read_left_posn == cur_read_left_posn)
		{
			rep_cnts[i_r] = rep_cnts[i_r-1] + 1;
		}
		else // This is a new fragment set its copy number to 1.
		{
			rep_cnts[i_r] = 1;
		}

		// Update prev. left posn.
		prev_read_left_posn = cur_read_left_posn;
    } // i_frag loop.

    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		if(rep_cnts[i_r] <= n_max_reps_per_posn)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		}
		else // This is a new fragment set its copy number to 1.
		{
			// Delete the pruned reads, otherwise they will be lost.
			delete_mapped_read(mapped_reads->at(i_r));

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Pruning %d repetition read @ %d, %c.\n", rep_cnts[i_r], mapped_reads->at(i_r)->base_index, mapped_reads->at(i_r)->strand);
			getc(stdin);
}
		}
    } // i_r loop.

	delete [] rep_cnts;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
	fprintf(stderr, "Pruned to %ld forward, %ld reverse strand fragments.\n", pruned_forward_reads->size(), pruned_reverse_reads->size());
}

bool sort_mapped_reads_per_5p(t_mapped_read* read1, t_mapped_read* read2)
{
	int frag1_5p = read_5p_accessor(read1);
	int frag2_5p = read_5p_accessor(read2);

	return(frag1_5p < frag2_5p);
}

int read_5p_accessor(void* obj_ptr)
{
	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;

	return(frag_obj_ptr->base_index);
}

//int read_3p_accessor(void* obj_ptr)
//{
//	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;
//
//	if(frag_obj_ptr->strand == 'F')
//	{
//		return(frag_obj_ptr->base_index+frag_obj_ptr->span-1);
//	}
//	else if(frag_obj_ptr->strand == 'R')
//	{
//		return(frag_obj_ptr->base_index);
//	}
//	else
//	{
//		fprintf(stderr, "The strand char for fragment object is %s @ %s(%d)\n", frag_obj_ptr->strand, __FILE__, __LINE__);
//		exit(0);
//		return(-1);
//	}
//}

