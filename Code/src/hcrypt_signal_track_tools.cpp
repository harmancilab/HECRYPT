#include "hcrypt_signal_track_tools.h"
#include "hcrypt_annot_region_tools.h"
#include "hcrypt_nomenclature.h"
#include "hcrypt_genomics_coords.h"
#include "hcrypt_xlog_math.h"
#include "hcrypt_utils.h"
#include "hcrypt_rng.h"
#include "hcrypt_ansi_string.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

#include "hcrypt_ansi_cli.h"
#include "hcrypt_config.h"

#include <ctype.h>
#include <string.h>
#include "hcrypt_histogram.h"
//#include "hcrypt_full_histogram.h"
#include "hcrypt_mapped_read_tools.h"
#include "hcrypt_ansi_string.h"
#include "hcrypt_seed_manager.h"
#include <algorithm>

bool __DUMP_SIGNAL_TRACK_MSGS__ = false;

bool sort_signal_nodes_per_increasing_signal(t_signal_node* node1, t_signal_node* node2)
{
	return(node1->signal < node2->signal);
}

double* invert_1based_signal(double* signal, int l_signal)
{
	double* inverted_signal = new double[l_signal + 2];
	memset(inverted_signal, 0, sizeof(double) * (l_signal + 2));

	for (int i = 1; i <= l_signal; i++)
	{
		inverted_signal[i] = signal[l_signal - i + 1];
	} // i loop.

	return(inverted_signal);
}

void reheader_signal_regions_BED(char* signal_regions_BED_fp, char* sample_ids_fp, char* op_fp)
{
	int n_loaded_samples = 0;
	vector<t_annot_region*>* all_sig_regs = load_signal_regs_BED(signal_regions_BED_fp, n_loaded_samples);
	fprintf(stderr, "Loaded %d regions and %d samples.\n", all_sig_regs->size(), n_loaded_samples);

	vector<char*>* sample_ids = buffer_file(sample_ids_fp);
	if (sample_ids->size() != n_loaded_samples)
	{
		fprintf(stderr, "Could not match the number of loaded samples and the matrix columns: %d, %d\n", sample_ids->size(), n_loaded_samples);
		exit(0);
	}

	fprintf(stderr, "Writing the reheadered file to %s\n", op_fp);

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_s = 0; i_s < sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < all_sig_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", all_sig_regs->at(i_reg)->chrom,
			translate_coord(all_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(all_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			all_sig_regs->at(i_reg)->name);

		double* cur_reg_sigs = (double*)(all_sig_regs->at(i_reg)->data);

		for (int i_s = 0; i_s < n_loaded_samples; i_s++)
		{
			// Make sure we do not count the first 4 columns.
			fprintf(f_op, "\t%lf", cur_reg_sigs[i_s]);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

void extract_signal_regions_per_sample_list_per_4th_col_signals(char* signal_reg_BED_fp, char* selected_sample_ids_list_fp, char* op_fp)
{
	int n_loaded_samples = 0;
	vector<t_annot_region*>* all_sig_regs = load_signal_regs_BED(signal_reg_BED_fp, n_loaded_samples);
	fprintf(stderr, "Loaded %d regions and %d samples.\n", all_sig_regs->size(), n_loaded_samples);

	char* header = load_header(signal_reg_BED_fp);
	vector<char*>* header_cols = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(header, "\t"));

	vector<char*>* selected_sample_ids = buffer_file(selected_sample_ids_list_fp);
	vector<char*>* config_sample_ids = new vector<char*>();
	vector<int>* per_config_sample_id_signal_reg_col_i = new vector<int>();
	for (int i_s = 0; i_s < selected_sample_ids->size(); i_s++)
	{
		config_sample_ids->push_back(t_string::copy_me_str(selected_sample_ids->at(i_s)));

		int cur_sample_i_col_i = t_string::get_i_str(header_cols, selected_sample_ids->at(i_s));
		if (cur_sample_i_col_i == header_cols->size())
		{
			fprintf(stderr, "Could not find %s\n", selected_sample_ids->at(i_s));
			exit(0);
		}

		per_config_sample_id_signal_reg_col_i->push_back(cur_sample_i_col_i);
	} // i_s loop.
	fprintf(stderr, "Extracting %d samples.\n", per_config_sample_id_signal_reg_col_i->size());

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_s = 0; i_s < per_config_sample_id_signal_reg_col_i->size(); i_s++)
	{
		fprintf(f_op, "\t%s", header_cols->at(per_config_sample_id_signal_reg_col_i->at(i_s)));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < all_sig_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", all_sig_regs->at(i_reg)->chrom,
			translate_coord(all_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(all_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			all_sig_regs->at(i_reg)->name);

		double* cur_reg_sigs = (double*)(all_sig_regs->at(i_reg)->data);

		for (int i_s = 0; i_s < per_config_sample_id_signal_reg_col_i->size(); i_s++)
		{
			// Make sure we do not count the first 4 columns.
			fprintf(f_op, "\t%lf", cur_reg_sigs[per_config_sample_id_signal_reg_col_i->at(i_s) - 4]);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

void concatenate_signal_regions(vector<char*>* signal_reg_BED_fps, char* op_fp)
{
	vector<t_annot_region*>** per_file_signal_regs = new vector<t_annot_region*>*[signal_reg_BED_fps->size() + 2];
	vector<int>* per_file_n_samples = new vector<int>();
	vector<char*>* per_file_headers = new vector<char*>();
	for (int i_f = 0; i_f < signal_reg_BED_fps->size(); i_f++)
	{
		int n_loaded_samples = 0;
		per_file_signal_regs[i_f] = load_signal_regs_BED(signal_reg_BED_fps->at(i_f), n_loaded_samples);
		fprintf(stderr, "%s: %d regions, %d samples.\n", signal_reg_BED_fps->at(i_f), per_file_signal_regs[i_f]->size(), n_loaded_samples);
		per_file_n_samples->push_back(n_loaded_samples);

		char* cur_file_header = load_header(signal_reg_BED_fps->at(i_f));
		per_file_headers->push_back(cur_file_header);
	} // i_f loop.

	vector<t_annot_region*>* agg_signal_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < per_file_signal_regs[0]->size(); i_reg++)
	{
		t_annot_region* cur_agg_reg = duplicate_region(per_file_signal_regs[0]->at(i_reg));
		t_annot_region** cur_agg_reg_per_file_regs = new t_annot_region*[signal_reg_BED_fps->size() + 2];
		memset(cur_agg_reg_per_file_regs, 0, sizeof(t_annot_region*) * (signal_reg_BED_fps->size() + 2));
		cur_agg_reg->data = cur_agg_reg_per_file_regs;
		agg_signal_regs->push_back(cur_agg_reg);
	} // i_reg loop.

	fprintf(stderr, "%d aggregate signal regions.\n", agg_signal_regs->size());

	for (int i_f = 0; i_f < signal_reg_BED_fps->size(); i_f++)
	{
		fprintf(stderr, "Adding %s\n", signal_reg_BED_fps->at(i_f));
		vector<t_annot_region*>* cur_intersects = intersect_annot_regions(agg_signal_regs, per_file_signal_regs[i_f], true);

		for (int i_int = 0; i_int < cur_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(cur_intersects->at(i_int)->data);
			t_annot_region* cur_agg_reg = int_info->src_reg;
			t_annot_region* cur_sig_reg = int_info->dest_reg;

			if (t_string::compare_strings(cur_agg_reg->name, cur_sig_reg->name))
			{
				t_annot_region** cur_agg_reg_per_file_regs = (t_annot_region**)(cur_agg_reg->data);
				cur_agg_reg_per_file_regs[i_f] = cur_sig_reg;
			}

			delete int_info;
		} // i_int loop.

		delete_annot_regions(cur_intersects);
	} // i_f loop

	fprintf(stderr, "Saving concatenated signal regions.\n");
	if (check_file(op_fp))
	{
		fprintf(stderr, "%s exists, will not overwrite.\n", op_fp);
		exit(0);
	}

	FILE* f_op = open_f(op_fp, "w");

	// Write the aggregate header.
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_f = 0; i_f < signal_reg_BED_fps->size(); i_f++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(per_file_headers->at(i_f), "\t");

		for (int i_t = 4; i_t < toks->size(); i_t++)
		{
			fprintf(f_op, "\t%s", toks->at(i_t)->str());
		} // i_t loop.

		t_string::clean_tokens(toks);
	} // i_f loop.
	fprintf(f_op, "\n");

	// Go over all the regions.
	for (int i_reg = 0; i_reg < agg_signal_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s", agg_signal_regs->at(i_reg)->chrom,
			translate_coord(agg_signal_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(agg_signal_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			agg_signal_regs->at(i_reg)->name);

		t_annot_region** cur_agg_reg_per_file_regs = (t_annot_region**)(agg_signal_regs->at(i_reg)->data);

		// Go over all the regions and dump.
		for (int i_f = 0; i_f < signal_reg_BED_fps->size(); i_f++)
		{
			if (cur_agg_reg_per_file_regs[i_f] == NULL)
			{
				fprintf(stderr, "Could not find sample %s for region %s:%d-%d\n", signal_reg_BED_fps->at(i_f),
					agg_signal_regs->at(i_reg)->chrom, agg_signal_regs->at(i_reg)->start, agg_signal_regs->at(i_reg)->end);
				exit(0);
			}

			double* cur_reg_signals = (double*)(cur_agg_reg_per_file_regs[i_f]->data);
			for (int i_s = 0; i_s < per_file_n_samples->at(i_f); i_s++)
			{
				fprintf(f_op, "\t%lf", cur_reg_signals[i_s]);
			} // i_s lo
		} // i_f loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}

// Use this to quantify total signal summaries over a range of samples.
// "Official" RPKM/TPM/RPM computation is in expression_tools and it is much more extendable.
void extract_signal_features_per_regions(char* argv[], int argc, char* config_fp, char* op_fp)
{
	t_config* config = new t_config(config_fp, "\t");
	t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");

	char feats_bed_fp[1000];
	if (!get_cli_value(config, cli, "feature_regs_fp", "-feature_regs_fp", feats_bed_fp, true))
	{
		fprintf(stderr, "Could not read: feature_regs_fp");
		exit(0);
	}
	
	// Read the element length norm. flag.
	char length_norm_flag_str[1000];
	if (!get_cli_value(config, cli, "length_normalization", "-length_normalization", length_norm_flag_str, true))
	{
		fprintf(stderr, "Could not read: length_normalization");
		exit(0);
	}
	bool length_norm_flag = (atoi(length_norm_flag_str) == 1);

	// Read the lib. size norm. flag.
	char lib_size_norm_flag_str[1000];
	if (!get_cli_value(config, cli, "lib_size_normalization", "-lib_size_normalization", lib_size_norm_flag_str, true))
	{
		fprintf(stderr, "Could not read: lib_size_normalization");
		exit(0);
	}
	bool lib_size_norm_flag = (atoi(lib_size_norm_flag_str) == 1);

	if (lib_size_norm_flag)
	{
		fprintf(stderr, "Performing lib. size normalization.\n");
	}
	else
	{
		fprintf(stderr, "NOT performing lib. size normalization.\n");
	}

	if (length_norm_flag)
	{
		fprintf(stderr, "Performing element length normalization.\n");
	}
	else
	{
		fprintf(stderr, "NOT performing element length normalization.\n");
	}

	// Load the set the signals.
	vector<t_annot_region*>* feat_regs_w_intervals = load_Regions_as_Interval(feats_bed_fp);
	for (int i_reg = 0; i_reg < feat_regs_w_intervals->size(); i_reg++)
	{
		double* cur_reg_per_sample_signals = new double[sample_entries->size() + 2];
		memset(cur_reg_per_sample_signals, 0, sizeof(double) * sample_entries->size());
		feat_regs_w_intervals->at(i_reg)->data = cur_reg_per_sample_signals;
	} // i_reg loop.

	vector<char*>* chr_ids = get_chr_ids(feat_regs_w_intervals);

	fprintf(stderr, "Extracting signal features on %d regions from %s.\n", feat_regs_w_intervals->size(), feats_bed_fp);

	// Total signal on each sample.
	double* per_sample_total_signal = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_signal_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal_in_regs, 0, sizeof(double) * sample_entries->size());
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_feat_regs = get_regions_per_chromosome(feat_regs_w_intervals, chr_ids->at(i_chr));
		fprintf(stderr, "Processing %d regions on %s\n", cur_chr_feat_regs->size(), chr_ids->at(i_chr));

		int signal_starting_col_i = 2;

		for (int i_s = 0; i_s < sample_entries->size(); i_s++)
		{
			// Pool the pileups/BGRs.
			int* cur_sample_cur_chr_covg = NULL;
			int l_pileup = -1;

			// Go over all the tracks for this sample and pool them.
			for (int col_i = signal_starting_col_i; col_i < sample_entries->at(i_s)->size(); col_i++)
			{
				int l_loaded_track_covg = 0;
				int* cur_sample_cur_chr_cur_track_covg = load_int_signal_covg_per_directory_chr_id(sample_entries->at(i_s)->at(col_i), chr_ids->at(i_chr), l_loaded_track_covg);

				// If the sample is not loaded, do not process.
				if (cur_sample_cur_chr_cur_track_covg == NULL)
				{
					fprintf(stderr, "Could not load signal from %s/%s; Skipping.\n", sample_entries->at(i_s)->at(col_i), chr_ids->at(i_chr));
					continue;
				}

				// Pool the current signal track.
				if (cur_sample_cur_chr_covg == NULL)
				{
					// Copy for the first track.
					l_pileup = l_loaded_track_covg + 1000;
					cur_sample_cur_chr_covg = new int[l_pileup];
					memset(cur_sample_cur_chr_covg, 0, sizeof(int) * l_pileup);
					for (int i = 1; i <= MIN(l_pileup, l_loaded_track_covg); i++)
					{
						cur_sample_cur_chr_covg[i] = cur_sample_cur_chr_cur_track_covg[i];
					} // i loop.
				}
				else
				{
					for (int i = 1; i <= MIN(l_pileup, l_loaded_track_covg); i++)
					{
						cur_sample_cur_chr_covg[i] += cur_sample_cur_chr_cur_track_covg[i];
					} // i loop.
				}

				delete[] cur_sample_cur_chr_cur_track_covg;
			} // col_i loop.

			// Update the total signal on this sample.
			for (int i = 1; i <= l_pileup; i++)
			{
				per_sample_total_signal[i_s] += cur_sample_cur_chr_covg[i];
			} // i loop.

			fprintf(stderr, "Sample %s: Total signal (chrom: %s): %.1f\n", sample_entries->at(i_s)->at(0), chr_ids->at(i_chr), per_sample_total_signal[i_s]);

			// Extract the total signal for all the regions on this chromosome.
			for (int i_reg = 0; i_reg < cur_chr_feat_regs->size(); i_reg++)
			{
				double* cur_reg_per_sample_signals = (double*)(cur_chr_feat_regs->at(i_reg)->data);

				vector<t_annot_region*>* cur_reg_int_regs = (vector<t_annot_region*>*)(cur_chr_feat_regs->at(i_reg)->intervals);
				vector<t_annot_region*>* merged_cur_reg_int_regs = merge_annot_regions(cur_reg_int_regs, 0);

				// Set the coverage to the dbl score of the region.
				double cur_reg_ints_coverage = coverage(merged_cur_reg_int_regs);
				cur_chr_feat_regs->at(i_reg)->dbl_score = cur_reg_ints_coverage;

				// Update the total signal in this region's intervals for the current sample; go over all the intervals of this element and update the signal.
				cur_reg_per_sample_signals[i_s] = 0;
				for (int i_int = 0; i_int < merged_cur_reg_int_regs->size(); i_int++)
				{
					for (int i = merged_cur_reg_int_regs->at(i_int)->start; i <= merged_cur_reg_int_regs->at(i_int)->end; i++)
					{
						if (i <= l_pileup)
						{
							cur_reg_per_sample_signals[i_s] += cur_sample_cur_chr_covg[i];
						}
						//else
						//{
						//	fprintf(stderr, "FATAL ERROR: Position further than chromosome length; %s: %d-%d (%d)\n",
						//		merged_cur_reg_int_regs->at(i_int)->chrom, merged_cur_reg_int_regs->at(i_int)->start, merged_cur_reg_int_regs->at(i_int)->end,
						//		l_pileup);

						//	exit(0);
						//}
					} // i loop.
				} // interval loop.

				  // Update the current total signal
				per_sample_total_signal_in_regs[i_s] += cur_reg_per_sample_signals[i_s];

				fprintf(stderr, "Sample %s: %s: %d exons (%d nucleotides in total): %.2f signal\r",
					sample_entries->at(i_s)->at(0),
					cur_chr_feat_regs->at(i_reg)->name,
					merged_cur_reg_int_regs->size(),
					(int)cur_reg_ints_coverage,
					cur_reg_per_sample_signals[i_s]);
				if (i_reg % 200 == 0)
				{
					fprintf(stderr, "\n");
				}

				delete_annot_regions(merged_cur_reg_int_regs);
			} // i_reg loop.

			delete[] cur_sample_cur_chr_covg;
		} // i_s loop.	
	} // i_chr loop

	  // Dump the normalized matrix of signal levels.
	FILE* f_op = open_f(op_fp, "w");

	// Dump the header.
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME ID TYPE\tFEATURE_LENGTH");

	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_op, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.

	fprintf(f_op, "\n");

	// Dump the signals for all the regions.
	for (int i_reg = 0; i_reg < feat_regs_w_intervals->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s\t%.1f", feat_regs_w_intervals->at(i_reg)->chrom,
			translate_coord(feat_regs_w_intervals->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(feat_regs_w_intervals->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			feat_regs_w_intervals->at(i_reg)->name,
			feat_regs_w_intervals->at(i_reg)->dbl_score);

		double* cur_reg_per_sample_signal = (double*)(feat_regs_w_intervals->at(i_reg)->data);
		for (int i_s = 0; i_s < sample_entries->size(); i_s++)
		{
			double cur_reg_normalized_signal = 0;

			if (lib_size_norm_flag)
			{
				cur_reg_normalized_signal = 1000 * 1000 * (cur_reg_per_sample_signal[i_s] / per_sample_total_signal[i_s]);
			}
			else
			{
				cur_reg_normalized_signal = cur_reg_per_sample_signal[i_s];
			}

			// Length normalization, if requested.
			if (length_norm_flag)
			{
				cur_reg_normalized_signal /= (feat_regs_w_intervals->at(i_reg)->dbl_score / 1000);
			}

			fprintf(f_op, "\t%lf", cur_reg_normalized_signal);
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);

	// Dump the signal statistics in the regions.
	FILE* f_per_sample_sig_stats = open_f("quant_signal_stats.txt", "w");
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_per_sample_sig_stats, "%.2f\t%.2f\t%.3f\n",
			per_sample_total_signal_in_regs[i_s], per_sample_total_signal[i_s],
			per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
	} // i_s loop.
	fclose(f_per_sample_sig_stats);
}

vector<t_annot_region*>* load_signal_regs_BED(char* signal_regions_BED_fp, int& n_loaded_samples, int matrix_starting_col_i)
{
	vector<t_annot_region*>* regs_w_lines = load_BED_with_line_information(signal_regions_BED_fp);
	vector<t_annot_region*>* signal_regs = new vector<t_annot_region*>();

	n_loaded_samples = -1;
	for (int i_reg = 0; i_reg < regs_w_lines->size(); i_reg++)
	{
		t_annot_region* cur_sig_reg = duplicate_region(regs_w_lines->at(i_reg));
		char* cur_reg_line = (char*)(regs_w_lines->at(i_reg)->data);

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");

		double* cur_sig_reg_sigs = new double[toks->size() + 2];
		for (int i_tok = matrix_starting_col_i; i_tok < toks->size(); i_tok++)
		{
			cur_sig_reg_sigs[i_tok - matrix_starting_col_i] = atof(toks->at(i_tok)->str());
		} // i_tok loop.
		cur_sig_reg->name = t_string::copy_me_str(toks->at(3)->str());
		cur_sig_reg->data = cur_sig_reg_sigs;

		if (n_loaded_samples == -1)
		{
			n_loaded_samples = toks->size() - matrix_starting_col_i;
		}
		else if (n_loaded_samples != toks->size() - matrix_starting_col_i)
		{
			fprintf(stderr, "Could not match the number of loaded samples: %d, %d\n", n_loaded_samples, toks->size() - matrix_starting_col_i);
			exit(0);
		}

		signal_regs->push_back(cur_sig_reg);

		t_string::clean_tokens(toks);
		delete[] cur_reg_line;
	} // i_reg loop.

	delete_annot_regions(regs_w_lines);

	return(signal_regs);
}

vector<t_annot_region*>* load_signal_regs_BED(char* signal_regions_BED_fp, int& n_loaded_samples)
{
	vector<t_annot_region*>* regs_w_lines = load_BED_with_line_information(signal_regions_BED_fp);
	vector<t_annot_region*>* signal_regs = new vector<t_annot_region*>();

	n_loaded_samples = -1;
	for (int i_reg = 0; i_reg < regs_w_lines->size(); i_reg++)
	{
		t_annot_region* cur_sig_reg = duplicate_region(regs_w_lines->at(i_reg));
		char* cur_reg_line = (char*)(regs_w_lines->at(i_reg)->data);

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");

		double* cur_sig_reg_sigs = new double[toks->size() + 2];
		for (int i_tok = 4; i_tok < toks->size(); i_tok++)
		{
			cur_sig_reg_sigs[i_tok - 4] = atof(toks->at(i_tok)->str());
		} // i_tok loop.
		cur_sig_reg->name = t_string::copy_me_str(toks->at(3)->str());
		cur_sig_reg->data = cur_sig_reg_sigs;

		if (n_loaded_samples == -1)
		{
			n_loaded_samples = toks->size() - 4;
		}
		else if (n_loaded_samples != toks->size() - 4)
		{
			fprintf(stderr, "Could not match the number of loaded samples: %d, %d\n", n_loaded_samples, toks->size() - 4);
			exit(0);
		}

		signal_regs->push_back(cur_sig_reg);

		t_string::clean_tokens(toks);
		delete[] cur_reg_line;
	} // i_reg loop.

	delete_annot_regions(regs_w_lines);

	return(signal_regs);
}

void remove_duplicate_samples_from_signal_matrix(char* signal_regs_bed_fp, char* op_signal_regs_bed_fp)
{
	int n_loaded_samples = 0;
	vector<t_annot_region*>* signal_regs = load_signal_regs_BED(signal_regs_bed_fp, n_loaded_samples);
	char* header = load_header(signal_regs_bed_fp);
	t_string_tokens* header_toks = t_string::tokenize_by_chars(header, "\t");
	vector<char*>* sample_ids = t_string::copy_tokens_2_strs(header_toks, 4, -1);
	if (sample_ids->size() != n_loaded_samples)
	{
		fprintf(stderr, "Loaded samples do not match the header.\n");
		exit(0);
	}

	vector<char*>* uniq_sample_ids = t_string::get_unique_entries(sample_ids);
	fprintf(stderr, "Mapping to %d unique samples out of %d samples.\n", uniq_sample_ids->size(), sample_ids->size());

	FILE* f_op = open_f(op_signal_regs_bed_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME ID TYPE");
	for (int i_s = 0; i_s < uniq_sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", uniq_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_reg = 0; i_reg < signal_regs->size(); i_reg++)
	{
		fprintf(f_op, "%s\t%d\t%d\t%s",
			signal_regs->at(i_reg)->chrom,
			translate_coord(signal_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(signal_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			signal_regs->at(i_reg)->name);

		double* cur_reg_signals = (double*)(signal_regs->at(i_reg)->data);

		for (int i_s = 0; i_s < uniq_sample_ids->size(); i_s++)
		{
			int cur_uniq_sample_i_s = t_string::get_i_str(sample_ids, uniq_sample_ids->at(i_s));

			if (cur_uniq_sample_i_s == sample_ids->size())
			{
				fprintf(stderr, "Could not find sample id %s\n", uniq_sample_ids->at(i_s));
				exit(0);
			}

			fprintf(f_op, "\t%lf", cur_reg_signals[cur_uniq_sample_i_s]);
		} // i_s loop.
		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);
}


void rank_normalize_signal_matrix(char* config_fp, char* signal_matrix_fp, char* op_fp)
{
	t_config* config = new t_config(config_fp, "\t");

	vector<char*>* matrix_lines = buffer_file(signal_matrix_fp);

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");

	vector<double>** per_sample_signals = new vector<double>*[sample_entries->size() + 2];
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		per_sample_signals[i_s] = new vector<double>();
	} // i_s loop.

	  // Parse the sample id's.
	t_string_tokens* sample_toks = t_string::tokenize_by_chars(matrix_lines->at(0), "\t");
	for (int i_t = 4; i_t < sample_toks->size(); i_t++)
	{
		fprintf(stderr, "Sample %d: %s\n", i_t - 4, sample_toks->at(i_t)->str());
	}

	for (int i_l = 1; i_l < matrix_lines->size(); i_l++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(matrix_lines->at(i_l), "\t");

		for (int i_t = 4; i_t < toks->size(); i_t++)
		{
			int i_s = i_t - 4;

			double cur_exp = atof(toks->at(i_t)->str());
			per_sample_signals[i_s]->push_back(cur_exp);
		} // i_t loop.
	} // i_l loop.

	fprintf(stderr, "Loaded %d regions with signals.\n", per_sample_signals[0]->size());

	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		vector<t_signal_node*>* nodes = new vector<t_signal_node*>();

		double cur_sample_total_sig = 0;
		for (int i_reg = 0; i_reg < per_sample_signals[i_s]->size(); i_reg++)
		{
			t_signal_node* node = new t_signal_node();
			node->signal = per_sample_signals[i_s]->at(i_reg);
			node->i_reg = i_reg;
			nodes->push_back(node);

			cur_sample_total_sig += per_sample_signals[i_s]->at(i_reg);
		} // i_reg loop.

		  //for (int i_reg = 0; i_reg < per_sample_signals[i_s]->size(); i_reg++)
		  //{
		  //	per_sample_signals[i_s]->at(i_reg) /= cur_sample_total_sig;
		  //	per_sample_signals[i_s]->at(i_reg) *= 1000;
		  //} // i_reg loop.

		sort(nodes->begin(), nodes->end(), sort_signal_nodes_per_increasing_signal);

		for (int rank_i = 0; rank_i < nodes->size(); rank_i++)
		{
			per_sample_signals[i_s]->at(nodes->at(rank_i)->i_reg) = rank_i;

			delete nodes->at(rank_i);
		} // rank_i loop.

		delete nodes;
	} // i_s loop.

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%s\n", matrix_lines->at(0));
	for (int i_l = 1; i_l < matrix_lines->size(); i_l++)
	{
		int i_reg = i_l - 1;

		t_string_tokens* toks = t_string::tokenize_by_chars(matrix_lines->at(i_l), "\t");
		fprintf(f_op, "%s\t%s\t%s\t%s", toks->at(0)->str(), toks->at(1)->str(), toks->at(2)->str(), toks->at(3)->str());
		for (int i_t = 4; i_t < toks->size(); i_t++)
		{
			int i_s = i_t - 4;

			fprintf(f_op, "\t%lf", per_sample_signals[i_s]->at(i_reg));
		} // i_t loop.

		fprintf(f_op, "\n");
	} // i_l loop.
	fclose(f_op);
}

double* load_signal_covg_per_signal_file(const char* cur_dat_fp,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded)
{
	double* covg_signal = NULL;

	reads_loaded = false;

	// Search for mapped reads.
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for BGR.
	if (check_file(cur_dat_fp) &&
		t_string::ends_with(cur_dat_fp, ".bgr.gz"))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	if (check_file(cur_dat_fp) &&
		t_string::ends_with(cur_dat_fp, ".bgr"))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

int* load_int_signal_covg_per_directory_chr_id(const char* dat_dir, char* chr_id, int& l_loaded_covg)
{
	int* covg_signal = NULL;

	// Search for pileup.
	char cur_dat_fp[1000];
	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for mapped reads.
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt", dat_dir, chr_id);
	if (!check_file(cur_dat_fp))
	{
		sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt.gz", dat_dir, chr_id);
	}

	int l_fragment = 0;
	if (check_file(cur_dat_fp))
	{
		//reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		double* dbl_covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			dbl_covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		int* covg_signal = new int[l_loaded_covg + 10];
		for (int i = 1; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.
		delete[] dbl_covg_signal;

		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/%s.bbgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bbgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bbgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		double* dbl_covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);

		covg_signal = new int[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (int)(dbl_covg_signal[i]);
		} // i loop.

		delete[] dbl_covg_signal;
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

double* load_signal_covg_per_directory_chr_id(const char* dat_dir,
	char* chr_id,
	int l_fragment,
	char* cur_dat_fp,
	int& l_loaded_covg,
	bool& reads_loaded)
{
	if (dat_dir == NULL)
	{
		return(NULL);
	}

	double* covg_signal = NULL;

	reads_loaded = false;

	// Search for mapped reads.
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for mapped reads.
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for pileup.
	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		int* int_covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		covg_signal = new double[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (double)(int_covg_signal[i]);
		} // i loop.

		delete[] int_covg_signal;

		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_allele_counts.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		int* int_covg_signal = load_coverage_per_compressed_pileup_file(cur_dat_fp, l_loaded_covg);
		covg_signal = new double[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = (double)(int_covg_signal[i]);
		} // i loop.

		delete[] int_covg_signal;

		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/%s.bbgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for bBGR.
	sprintf(cur_dat_fp, "%s/%s.bbgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_bBGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s.bin", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_binary_profile(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_binary_profile(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_logR.bin.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_binary_profile(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

double* concat_profile2_2_profile1(double* prof1, int l_prof1, 
									double* prof2, int l_prof2, 
									int& l_new_track)
{
	double* conc_prof = NULL;
	if(prof2 == NULL)
	{
		l_new_track = l_prof1;

		conc_prof = new double[l_new_track + 1];

		for(int i = 1; i <= l_prof1; i++)
		{
			conc_prof[i] = prof1[i];
		} // i loop.
	}
	else if(prof1 == NULL)
	{
		l_new_track = l_prof2;

		conc_prof = new double[l_new_track + 1];

		for(int i = 1; i <= l_prof2; i++)
		{
			conc_prof[i] = prof2[i];
		} // i loop.

	}
	else
	{
		l_new_track = l_prof1 + l_prof2;

		conc_prof = new double[l_new_track + 1];
		for(int i = 1; i <= l_prof1; i++)
		{
			conc_prof[i] = prof1[i];
		} // i loop.

		for(int i = 1; i <= l_prof2; i++)
		{
			conc_prof[i + l_prof1] = prof2[i];
		} // i loop.
	}

	return(conc_prof);
}

void dump_per_nucleotide_uchar_binary_profile(unsigned char* signal_profile, int l_profile, char* op_fp)
{
	// Dump the per nucleotide signal profile.
	FILE* f_op = open_f(op_fp, "wb");

	// Write data length to first couple bytes.
	fwrite(&(l_profile), sizeof(int), 1, f_op);

	// Dump the data: Dump 0 based data.
	fwrite(&(signal_profile[1]), sizeof(unsigned char), l_profile+1, f_op);

	fclose(f_op);
}

unsigned char* load_per_nucleotide_binary_uchar_profile(char* binary_per_nucleotide_profile_fp, int& l_profile)
{
	FILE* f_prof = open_f(binary_per_nucleotide_profile_fp, "rb");

	// Read the profile length.
	int l_data = 0;
	fread(&l_data, sizeof(int), 1, f_prof);
	l_profile = l_data;

	// Read the data.
	unsigned char* signal_profile_buffer = new unsigned char[l_profile+2];

if(__DUMP_SIGNAL_TRACK_MSGS__)
	fprintf(stderr, "Loading %d data values.\n", l_profile);

	// Following is to use the codebase indexing: 1 based indexing.
	fread(&(signal_profile_buffer[1]), sizeof(char), l_profile+1, f_prof);
	
	fclose(f_prof);

	return(signal_profile_buffer);
}

double* get_block_permute_profile(double* profile_data, t_rng* rng, int l_profile, int& l_permuted, int l_block)
{
	//int l_cur_profile;
	//double* cur_profile_data = load_per_nucleotide_binary_profile(profile_fp, l_cur_profile);

	//  Count the blocks.
	int n_blocks = (int)(l_profile / l_block);
	fprintf(stderr, "There are %d blocks in total.\n", n_blocks);

	// Permute the block ordering.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me());

	//vector<int>* permuted_block_i = rng->permute_indices(n_blocks, n_blocks);
	vector<int>* permuted_block_i = rng->fast_permute_indices(0, n_blocks);

	double* permuted_profile_data = new double[l_profile + 2];
	l_permuted = 1;
	for(int block_i = 0; block_i < (int)permuted_block_i->size(); block_i++)
	{
		int cur_block_start = (permuted_block_i->at(block_i) * l_block)+1;
		int cur_block_end = (l_profile > cur_block_start+l_block)?(cur_block_start+l_block):(l_profile);

		for(int i = cur_block_start; i < cur_block_end; i++)
		{
			permuted_profile_data[l_permuted] = profile_data[i];
			l_permuted++;
		} // i loop.
	} // block_i loop.

	// Dump the profile.
	return(permuted_profile_data);
}

void get_profile_extrema(double* profile, int l_profile, double& prof_min, double& prof_max)
{
	prof_min = 1000*1000;
	prof_max = -1000*1000;
	for(int i = 1; i <= l_profile; i++)
	{
		if(profile[i] > prof_max)
		{
			prof_max = profile[i];
		}

		if(profile[i] < prof_min)
		{
			prof_min = profile[i];
		}
	} // i loop.
}

double* get_zero_indexed_per_one_indexed_data(double* one_indexed_data, int l_profile)
{
	double* zero_indexed_data = new double[l_profile];
	for(int i = 1; i <= l_profile; i++)
	{
		zero_indexed_data[i-1] = one_indexed_data[i];
	} // i loop.

	return(zero_indexed_data);
}

double* get_one_indexed_per_zero_indexed_data(double* zero_indexed_data, int l_profile)
{
	double* one_indexed_data = new double[l_profile+2];
	for(int i = 0; i < l_profile; i++)
	{
		one_indexed_data[i+1] = zero_indexed_data[i];
	} // i loop.

	return(one_indexed_data);
}

double* quantize_per_nucleotide_profiles(double* profile, int l_profile, vector<double>* thresholds, vector<double>* quantized_vals)
{
	//int l_profile = 0;
	//double* profile = load_per_nucleotide_binary_profile(prof_fp, l_profile);

	//vector<char*>* thresh_val_lines = buffer_file(thresh_val_fp);

	//vector<double>* quantized_vals = new vector<double>();
	//vector<double>* thresholds = new vector<double>();
	//for(int i_l = 0; i_l < thresh_val_lines->size(); i_l++)
	//{
	//	double cur_quantized_val;
	//	double cur_thresh;
	//	if(sscanf(thresh_val_lines->at(i_l), "%lf %lf", &cur_thresh, &cur_quantized_val) != 2)
	//	{
	//		fprintf(stderr, "Could not parse the line: %s\n", thresh_val_lines->at(i_l));
	//		exit(0);
	//	}

	//	quantized_vals->push_back(cur_quantized_val);
	//	thresholds->push_back(cur_thresh);
	//	fprintf(stderr, "%lf -> %lf\n", cur_thresh, cur_quantized_val);
	//} // i_l loop.

	if(thresholds->size() != quantized_vals->size())
	{
		fprintf(stderr, "The size of thresholds is not the same as size of quantized values.\n");
		exit(0);
	}

	double* quantized_profile = new double[l_profile+2];
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_prof_val = profile[i];

		// Find the largest threshold that this value is greater than or equal to.
		bool quantized_current_val = false;
		for(int i_th = thresholds->size()-1;
			!quantized_current_val && i_th >= 0; 
			i_th--)
		{
			if(cur_prof_val >= thresholds->at(i_th))
			{
				quantized_current_val = true;
				quantized_profile[i] = quantized_vals->at(i_th);
			}
		} // i_th loop.

		if(!quantized_current_val)
		{
			fprintf(stderr, "Could not quantize the current value: %lf\n", profile[i]);
			exit(0);
		}
	} // i loop.

	return(quantized_profile);
}

double* copy_profile(double* signal_profile, int l_profile)
{
	double* copy_prof = new double[l_profile+2];

	for(int i_sig = 1; i_sig <= l_profile; i_sig++)
	{
		copy_prof[i_sig] = signal_profile[i_sig];
	} // i_sig loop.

	return(copy_prof);
}

double* extract_one_indexed_profile_per_profile_stranded(double* signal_profile_buffer, int l_profile, int start, int end, char strand, int& l_extracted_profile)
{
	double* extracted_prof = new double[end - start + 3];
	for (int i_sig = start; i_sig <= end; i_sig++)
	{
		extracted_prof[i_sig - start] = 0;
	} // i_sig loop.

	// Count the extracted profile.
	l_extracted_profile = 1;

	if (strand == '+' ||
		strand == 'F')
	{
		for (int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
		{
			extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
			l_extracted_profile++;
		} // i_sig loop.
	}
	else
	{
		//for (int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
		for (int i_sig = MIN(l_profile, end); i_sig >= start; i_sig--)
		{
			extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
			l_extracted_profile++;
		} // i_sig loop.
	}
	// Profiles are one based.
	l_extracted_profile--;

	return(extracted_prof);
}


double* extract_one_indexed_profile_per_profile(double* signal_profile_buffer, int l_profile, int start, int end, int& l_extracted_profile)
{
	double* extracted_prof = new double[end - start + 3];	
	for(int i_sig = start; i_sig <= end; i_sig++)
	{
		extracted_prof[i_sig - start] = 0;
	} // i_sig loop.
	
	// Count the extracted profile.
	l_extracted_profile = 1;
	for(int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
	{
		extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
		l_extracted_profile++;
	} // i_sig loop.

	// Profiles are one based.
	l_extracted_profile--;

	return(extracted_prof);
}

// The signal profile is 1 based, consistent with the codebase indexing.
void dump_bedGraph_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* op_fp)
{
	FILE* f_op = NULL;
	bool concatting = false;
	if(check_file(op_fp))
	{
		fprintf(stderr, "%s exists, concatting.\n", op_fp);
		f_op = open_f(op_fp, "a");
		concatting = true;
	}
	else
	{
		f_op = open_f(op_fp, "w");
	}

	// Get the bedgraph for the current profile.
	int i_nuc = 1; 
	double cur_height = signal_profile_buffer[i_nuc];
	int cur_block_start_i = i_nuc;
	i_nuc++;
	while(1)
	{
		// Find the point where the height changes: The end goes till it's equal to the profile length since the profile is 1-based.
		while(i_nuc <= l_profile)
		{
			// Wait till there is a change in the height, which marks the start of a new block.
			if(cur_height != signal_profile_buffer[i_nuc])
			{
				break;
			}

			i_nuc++;
		} // i_nuc loop.

		// At this point, either this is the end of the profile, or there was a change in the height, either way, this was the end of the current block. Definitely dump it.
		if(cur_height != signal_profile_buffer[i_nuc])
		{
			// Dump the current block.
			fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom, 
				translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(i_nuc-1, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
				cur_height);

			// Update the new height and new start.
			cur_height = signal_profile_buffer[i_nuc];
			
			// Current position starts the next block.
			cur_block_start_i = i_nuc; 
		}

		// If the above block end was the end of the whole profile, we are done, otherwise continue to the next block.
		if(i_nuc > l_profile)
		{
			break;
		}

		//i_nuc++;
	} // i_nuc loop.

	//fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom, 
	//			translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//			translate_coord(i_nuc-1, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
	//			cur_height);

	close_f(f_op, op_fp);
}

// The signal profile is 1 based, consistent with the codebase indexing: This is the latest version of bBGR dumping.
void dump_bBGR_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* bbgr_op_fp)
{
	FILE* f_bbgr_op = open_f(bbgr_op_fp, "wb");

	// Get the bedgraph for the current profile.
	int i_nuc = 1;
	double cur_height = signal_profile_buffer[i_nuc];
	int cur_block_start_i = i_nuc;
	i_nuc++;
	while (1)
	{
		// Find the point where the height changes: The end goes till it's equal to the profile length since the profile is 1-based.
		while (i_nuc <= l_profile)
		{
			// Wait till there is a change in the height, which marks the start of a new block.
			if (cur_height != signal_profile_buffer[i_nuc])
			{
				break;
			}

			i_nuc++;
		} // i_nuc loop.

		  // At this point, either this is the end of the profile, or there was a change in the height, either way, this was the end of the current block. Definitely dump it.
		if (cur_height != signal_profile_buffer[i_nuc])
		{
			int start = translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			int end = translate_coord(i_nuc - 1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			//// Dump the current block.
			//fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom,
			//	translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			//	translate_coord(i_nuc - 1, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			//	cur_height);

			// Dump the current block.
			fwrite(&start, sizeof(int), 1, f_bbgr_op);
			fwrite(&end, sizeof(int), 1, f_bbgr_op);
			fwrite(&cur_height, sizeof(double), 1, f_bbgr_op);

			// Update the new height and new start.
			cur_height = signal_profile_buffer[i_nuc];

			// Current position starts the next block.
			cur_block_start_i = i_nuc;
		}

		// If the above block end was the end of the whole profile, we are done, otherwise continue to the next block.
		if (i_nuc > l_profile)
		{
			break;
		}

		//i_nuc++;
	} // i_nuc loop.

	fclose(f_bbgr_op);

	char comp_bbgr_op_fp[1000];
	sprintf(comp_bbgr_op_fp, "%s.gz", bbgr_op_fp);
	compressFile(bbgr_op_fp, comp_bbgr_op_fp);
	delete_file(bbgr_op_fp);
}

void dump_bedGraph_per_bBGR_v1(char* bbgr_fp, char* chrom, char* op_fp)
{
	int l_profile = 0;
	double* profile = load_per_nucleotide_bBGR_v1_track(bbgr_fp, l_profile);

	dump_bedGraph_per_per_nucleotide_binary_profile(profile, l_profile, chrom, op_fp);
}

void dump_bedGraph_per_bBGR(char* bbgr_fp, char* chrom, char* op_fp)
{
	int l_profile = 0;
	double* profile = load_per_nucleotide_bBGR_track(bbgr_fp, l_profile);

	dump_bedGraph_per_per_nucleotide_binary_profile(profile, l_profile, chrom, op_fp);
}

void dump_bBGR_per_per_bedGraph(char* bgr_fp, char* bbgr_op_fp)
{
	FILE* f_bgr = open_f(bgr_fp, "r");

	FILE* f_bbgr_op = open_f(bbgr_op_fp, "wb");
	while (1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if (cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if (sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		//fwrite(cur_chrom, 1, 20, f_bbgr_op);
		fwrite(&start, sizeof(int), 1, f_bbgr_op);
		fwrite(&end, sizeof(int), 1, f_bbgr_op);
		fwrite(&cur_sig, sizeof(double), 1, f_bbgr_op);
	} // bgr reading loop.

	close_f(f_bbgr_op, bbgr_op_fp);
	close_f(f_bgr, bgr_fp);
}

double* load_per_nucleotide_bBGR_track(char* bgr_fp, int& l_profile)
{
	FILE* f_bbgr = open_f(bgr_fp, "rb");

	// Initialize the signal profile.
	int l_max = 300 * 1000 * 1000;
	double* signal_profile = new double[l_max + 1];
	memset(signal_profile, 0, sizeof(double) * l_max);

	// Go over all the input file and process all the lines.
	while (1)
	{
		int i_chr = 0;
		int start = 0;
		int end = 0;
		double cur_sig = 0;
		if (fread(&start, sizeof(int), 1, f_bbgr) == 0)
		{
			break;
		}

		fread(&end, sizeof(int), 1, f_bbgr);
		fread(&cur_sig, sizeof(double), 1, f_bbgr);

		l_profile = end + 10;

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for (int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if (i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d\n", i_nuc);
				exit(0);
			}
			signal_profile[i_nuc] += cur_sig;
		} // i_nuc loop.
	} // file reading loop.	

	// Close the bedgraph file.
	close_f(f_bbgr, bgr_fp);

	return(signal_profile);
}

double* load_per_nucleotide_bBGR_v1_track(char* bgr_fp, int& l_profile)
{
	FILE* f_bbgr = open_f(bgr_fp, "rb");

	// Initialize the signal profile.
	int l_max = 300 * 1000 * 1000;
	double* signal_profile = new double[l_max + 1];
	memset(signal_profile, 0, sizeof(double) * l_max);

	// Go over all the input file and process all the lines.
	while (1)
	{
		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if (fread(cur_chrom, 1, 20, f_bbgr) == 0)
		{
			break;
		}

		fread(&start, sizeof(int), 1, f_bbgr);
		fread(&end, sizeof(int), 1, f_bbgr);
		fread(&cur_sig, sizeof(double), 1, f_bbgr);

		l_profile = end + 10;

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for (int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if (i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d\n", i_nuc);
				exit(0);
			}
			signal_profile[i_nuc] = cur_sig;
		} // i_nuc loop.
	} // file reading loop.	

	  // Close the bedgraph file.
	close_f(f_bbgr, bgr_fp);

	return(signal_profile);
}

double* load_per_nucleotide_BGR_track(const char* bgr_fp, int& l_profile)
{
	FILE* f_bgr = open_f(bgr_fp, "r");

	// Initialize the signal profile.
	int l_max = 300 * 1000 * 1000;
	double* signal_profile = new double[l_max + 1];
	memset(signal_profile, 0, sizeof(double) * l_max);

	// Go over all the input file and process all the lines.
	while (1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if (cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if (sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		l_profile = end + 10;

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for (int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if (i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d (%s)\n", i_nuc, cur_bgr_line);
				exit(0);
			}

			// Note that this pools if there are overlapping positions in the bedgraph file.
			signal_profile[i_nuc] += cur_sig;
		} // i_nuc loop.

		delete[] cur_bgr_line;
	} // file reading loop.	

	// Close the bedgraph file.
	close_f(f_bgr, bgr_fp);

	return(signal_profile);
}

void dump_per_nucleotide_binary_profile_per_bedgraph(const char* bgr_fp, bool dump_binary, const char* op_fp)
{
	FILE* f_bgr = open_f(bgr_fp, "r");

	// Initialize the signal profile.
	int l_max = 300*1000*1000;
	double* signal_profile = new double[l_max+1];
	for(int i_nuc = 0; i_nuc <= l_max; i_nuc++)
	{
		signal_profile[i_nuc] = 0.0;
	} // i_nuc loop.

	// Go over all the input file and process all the lines.
	fprintf(stderr, "Dumping the profile to %s.\n", op_fp);
	while(1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if(cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if(sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for(int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if(i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d (%s)\n", i_nuc, cur_bgr_line);
				exit(0);
			}
			signal_profile[i_nuc] = cur_sig;
		} // i_nuc loop.

		delete [] cur_bgr_line;
	} // file reading loop.	

	// Close the bedgraph file.
	close_f(f_bgr, bgr_fp);

	// Get the end of the signal.
	int l_data = l_max;
	while(signal_profile[l_data] == 0.0)
	{
		l_data--;
	} // i_nuc loop.

	fprintf(stderr, "Signal length is %d, dumping the per nucleotide profile.\n", l_data);

	if(dump_binary)
	{
		dump_per_nucleotide_binary_profile(signal_profile, l_data, op_fp);
		return;
	}

	// Dump the per nucleotide signal profile.
	FILE* f_op = NULL;
	
	fprintf(stderr, "Saving Plain + Compressing.\n");
	f_op = open_f(op_fp, "w");
	
	fprintf(f_op, "%d\n",  l_data);
	
	for(int i_nuc = 1; i_nuc <= l_data; i_nuc++)
	{
		if(dump_binary)
		{
			fwrite(&(signal_profile[i_nuc]), sizeof(double), 1, f_op);
		}
		else
		{
			fprintf(f_op, "%lf ",  signal_profile[i_nuc]);
		}
	} // i_nuc loop.

	close_f(f_op, op_fp);
}

void dump_per_nucleotide_binary_profile(double* signal_profile, int l_profile, const char* op_fp)
{
	// Dump the per nucleotide signal profile.
	FILE* f_op = open_f(op_fp, "wb");

	// Write data length to first couple bytes.
	fwrite(&(l_profile), sizeof(int), 1, f_op);

	// Dump the data: Dump 0 based data.
	fwrite(&(signal_profile[1]), sizeof(double), l_profile, f_op);

	close_f(f_op, op_fp);
}

double* load_per_nucleotide_binary_profile(char* binary_per_nucleotide_profile_fp, int& l_profile)
{
	FILE* f_prof = open_f(binary_per_nucleotide_profile_fp, "rb");

	if (f_prof == NULL)
	{
		fprintf(stderr, "Could not open %s, unexpected extension other than bin/bin.gz?\n", binary_per_nucleotide_profile_fp);
		exit(0);
	}

	// Read the profile length.
	int l_data = 0;
	fread(&l_data, sizeof(int), 1, f_prof);
	l_profile = l_data;

	// Read the data.
	double* signal_profile_buffer = new double[l_profile+2];
	fprintf(stderr, "Loading %d data values.\n", l_profile);

	// Following is to use the codebase indexing: 1 based indexing.
	fread(&(signal_profile_buffer[1]), sizeof(double), l_profile, f_prof);
	
	close_f(f_prof, binary_per_nucleotide_profile_fp);

	return(signal_profile_buffer);
}

void exclude_regions_from_signal_profiles(double* signal_profile, int l_profile, vector<t_annot_region*>* regions_2_exclude, double* pruned_signal_profile, int& l_pruned_profile)
{
	sort(regions_2_exclude->begin(), regions_2_exclude->end(), sort_regions);

	bool* exclusion_profile = new bool[l_profile+2];
	for(int i = 1; i <= l_profile; i++)
	{
		exclusion_profile[i] = false;
	} // i loop.

	for(int i_reg = 0; i_reg < (int)regions_2_exclude->size(); i_reg++)
	{
		for(int i = regions_2_exclude->at(i_reg)->start; i <= regions_2_exclude->at(i_reg)->end; i++)
		{
			if(i <= l_profile)
			{
				exclusion_profile[i] = true;
			}
		} // i loop.
	} // i_reg loop.

	l_pruned_profile = 1;

	for(int i = 1; i <= l_profile; i++)
	{
		if(exclusion_profile[i])
		{
		}
		else
		{
			pruned_signal_profile[l_pruned_profile] = signal_profile[i];
			l_pruned_profile++;
		}
	} // i loop.

	// Must move the pruned profile index back by one.
	l_pruned_profile--;

	fprintf(stderr, "Pruned the signal of %d values to %d values.\n", l_profile, l_pruned_profile);
	delete [] exclusion_profile;
}

void get_log_plus_one_profile(double* signal_profile, double base, int l_profile)
{
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_log_value = xlog(signal_profile[i]+1.0)/xlog(base);
		signal_profile[i] = cur_log_value;
	} // i loop.
}

void floorize_profile(double* signal_profile, int l_profile)
{
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_floor_val = floor(signal_profile[i]);
		signal_profile[i] = cur_floor_val;
	} // i loop.
}

double* get_zero_profile(int l_profile)
{
	double* cur_profile = new double[l_profile+2];
	for(int i = 0; i <= l_profile; i++)
	{
		cur_profile[i] = 0.0;
	} // i loop.

	return(cur_profile);
}

