#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"
#include "time.h"
#include <chrono>
#include <sys/resource.h>   // check the memory
#include <random>
#include <iomanip>          // std::setprecision

#include "../../external/UTMSR/src/param.h"
#include "../../external/UTMSR/src/utils.h"
#include "../../external/UTMSR/src/HEmpute_bfv.h"
#include "../../external/UTMSR/src/HEmpute_ckks.h"
#include "../../external/UTMSR/src/TestHEmpute.h"

#include "../../external/UTMSR/src/utils_data.h"
#include "../../external/UTMSR/src/thread.h"

#include "hcrypt_variation_tools.h"
#include "hcrypt_annot_region_tools.h"
#include "hcrypt_gff_utils.h"
#include "hcrypt_histogram.h"
#include "hcrypt_genomics_coords.h"
#include "hcrypt_utils.h"
#include "hcrypt_ansi_string.h"
#include "hcrypt_signal_track_tools.h"

#include <vector>
#include "seal/seal.h"

using namespace std;
using namespace seal;

#include "hcrypt_hecrypt_utils.h"

#define DEBUG false
#define fast true   // using the fast multiplications between single-precision numbers

// This function checks wheter the tag at the position is the closest leftmost tag, including target overlapping tags.
bool hempute_OLS_closest_left_tag_checker(vector<t_annot_region*>* tag_var_regs, int tag_var_i,
	int target_start_posn)
{
	if (tag_var_regs->at(tag_var_i)->start <= target_start_posn &&
		tag_var_regs->at(tag_var_i + 1)->start > target_start_posn)
	{
		return(true);
	}

	return(false);
}

// This function checks wheter the tag at the position is the closest leftmost tag, wrapper for regions as arguments.
bool hempute_OLS_closest_left_tag_checker(vector<t_annot_region*>* tag_var_regs, int tag_var_i,
	vector<t_annot_region*>* target_var_regs, int target_var_i)
{
	return(hempute_OLS_closest_left_tag_checker(tag_var_regs, tag_var_i,
		target_var_regs->at(target_var_i)->start));
}

void get_R2_per_imputed_genotypes(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp, char* stats_op_fp)
{
	vector<t_annot_region*>* known_genotype_regs = load_variant_signal_regions_wrapper(known_genotypes_fp, known_sample_ids_list_fp);
	vector<char*>* known_sample_ids = buffer_file(known_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", known_genotype_regs->size(), known_sample_ids->size());
	t_restr_annot_region_list* restr_known_genotype_regs = restructure_annot_regions(known_genotype_regs);

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", imputed_sample_ids->size());

	int n_loaded_samples = 0;
	vector<t_annot_region*>* imputed_genotype_regs = load_signal_regs_BED(imputed_genotypes_fp, n_loaded_samples);
	fprintf(stderr, "Loaded %d imputed variants.\n", imputed_genotype_regs->size());
	for (int i_reg = 0; i_reg < imputed_genotype_regs->size(); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
		double* cur_reg_imp_sig = (double*)(imputed_genotype_regs->at(i_reg)->data);
		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = cur_reg_imp_sig;
		imputed_genotype_regs->at(i_reg)->data = cur_reg_info;
	} // i_reg loop.
	t_restr_annot_region_list* restr_imputed_genotype_regs = restructure_annot_regions(imputed_genotype_regs);

	// Set the mapping between sample id's.
	vector<int>* imp_2_known_sample_i = new vector<int>();
	int n_matched_samples = 0;
	for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
	{
		int cur_imp_known_i = t_string::get_i_str(known_sample_ids, imputed_sample_ids->at(imp_i));
		if (cur_imp_known_i < known_sample_ids->size())
		{
			n_matched_samples++;
		}
		imp_2_known_sample_i->push_back(cur_imp_known_i);
	} // imp_i loop.
	fprintf(stderr, "Matched %d samples.\n", n_matched_samples);

	FILE* f_op = open_f(stats_op_fp, "w");
	for (int i_chr = 0; i_chr < restr_known_genotype_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Intersecting known regions with imputed regions.\n");
		
		// Intersect and process.
		double** known_imp_sample_geno = new double*[2];
		known_imp_sample_geno[0] = new double[n_matched_samples + 2];
		known_imp_sample_geno[1] = new double[n_matched_samples + 2];
		vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);
		fprintf(stderr, "Found %d intersections\n", intersects->size());
		int n_processed_imputed_targets = 0;
		for (int i_int = 0; i_int < intersects->size(); i_int++)
		{
			if (i_int % 1000 == 0)
			{
				fprintf(stderr, "@ %d. intersect         \r", i_int);
			}

			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* imp_reg = int_info->src_reg;
			t_annot_region* known_reg = int_info->dest_reg;

			void** known_reg_info = (void**)(known_reg->data);
			char* known_reg_geno = (char*)(known_reg_info[0]);

			void** imp_reg_info = (void**)(imp_reg->data);
			double* imp_reg_geno = (double*)(imp_reg_info[0]);

			// We do not enforce region name matching.
			//if (t_string::compare_strings(imp_reg->name, known_reg->name) &&
			if(imp_reg->score == 0)
			{
				imp_reg->score = 1;
				n_processed_imputed_targets++;

				double total_n_non_refs = 0;
				double n_matching_non_refs = 0;
				double n_matching_all = 0;
				double n_all = 0;

				int geno_i = 0;
				for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
				{
					if (imp_2_known_sample_i->at(imp_i) < known_sample_ids->size())
					{
						double known_geno = (double)(known_reg_geno[imp_2_known_sample_i->at(imp_i)]);
						double imp_geno = (double)(imp_reg_geno[imp_i]);
						known_imp_sample_geno[0][geno_i] = known_geno;
						known_imp_sample_geno[1][geno_i] = imp_geno;
						geno_i++;

						// Update non-ref concordance.
						if (known_geno > 0)
						{
							if (known_geno == imp_geno)
							{
								n_matching_non_refs++;
							}

							total_n_non_refs++;
						}

						// Update all concordance.
						if (known_geno == imp_geno)
						{
							n_matching_all++;
						}
						n_all++;
					} // matching check.
				} // imp_i loop.			

				double cur_geno_corr = 0;
				get_correlation(known_imp_sample_geno[0], known_imp_sample_geno[1], geno_i, cur_geno_corr);

				fprintf(f_op, "%s\t%d\t%d\t%s\t%.4f\t%.0f\t%.0f\t%.0f\t%.0f\n", imp_reg->chrom,
					translate_coord(imp_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(imp_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					imp_reg->name, cur_geno_corr * cur_geno_corr,
					n_matching_all, n_all,
					n_matching_non_refs, total_n_non_refs);
			} // overlapping region name comparison.
		} // i_int loop.	
	} // i_chr loop.
	fclose(f_op);
	fprintf(stderr, "\nDone.\n");
} // get_R2_per_imputed_genotypes function.

void compare_LMSE_vs_secure_eval(char* lmse_sig_fp, char* secure_val_sig_fp, char* sample_ids_list_fp)
{
	int n_lmse_samples = 0;
	vector<t_annot_region*>* lmse_regs = load_signal_regs_BED(lmse_sig_fp, n_lmse_samples);

	int n_secure_eval_samples = 0;
	vector<t_annot_region*>* secure_eval_regs = load_signal_regs_BED(secure_val_sig_fp, n_secure_eval_samples);

	if (n_lmse_samples != n_secure_eval_samples)
	{
		fprintf(stderr, "Sample sizes not matchng, %d/%d\n", n_lmse_samples, n_secure_eval_samples);
		exit(0);
	}

	fprintf(stderr, "Loaded %d LMSE and %d secure eval regions.\n", lmse_regs->size(), secure_eval_regs->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(lmse_regs, secure_eval_regs, false);
	fprintf(stderr, "Found %d intersects.\n", intersects->size());

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* lmse_reg = int_info->src_reg;
		t_annot_region* secure_eval_reg = int_info->dest_reg;

		double* lmse_sig = (double*)(lmse_reg->data);
		double* secure_eval_sig = (double*)(secure_eval_reg->data);

		for (int i_s = 0; i_s < n_lmse_samples; i_s++)
		{
			if (fabs((lmse_sig[i_s] * 2) - secure_eval_sig[i_s]) > 0.1)
			{
				fprintf(stderr, "%d: %d: %.4f/%.4f\n", 
						intersects->at(i_int)->start, 
						i_s, 
						lmse_sig[i_s] / 50, secure_eval_sig[i_s]);
			}
		} // i_s loop.
	} // i_int loop.
}

// This is a testing method for decryption: This decrypts the testing genotypes. Need a new data for target variants.
void generic_decrypt_genotypes(char* encrypted_geno_matrix_fp,
										char* tag_geno_coords_fp,
										char* sample_ids_list_fp,
										char* key_prefix,
										const char* decrypted_geno_op_fp)
{
	fprintf(stderr, "Decrypting %s\n", encrypted_geno_matrix_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	vector<t_annot_region*>* tag_variant_coord_regs = load_BED(tag_geno_coords_fp);
	fprintf(stderr, "Loaded %d tag variant coordinates.\n", tag_variant_coord_regs->size());

	EncryptionParameters parms(scheme_type::BFV);
	size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
	int logq = 27;

	// The plain_modulus does not play much of a role;
	// we choose some reasonable value t > results so that it satisfies (results mod t) = results
	size_t plaintext_modulus = (1 << 10);

	double Xscale = 2.0;
	double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
	double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg

	parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { logq }));
	parms.set_plain_modulus(plaintext_modulus);

	auto context = SEALContext::Create(parms);
	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	// Write the public/private keys.
	fprintf(stderr, "Loading keys.\n");
	char key_fp[1000];
	sprintf(key_fp, "%s.public_key", key_prefix);

	if (check_file(key_fp))
	{
		fprintf(stderr, "Loading %s\n", key_fp);
		std::ifstream pub_key_file;
		pub_key_file.open(key_fp, std::ifstream::in);
		public_key.load(context, pub_key_file);
		pub_key_file.close();
	}
	else
	{
		fprintf(stderr, "Public key does not exist @ %s -- Make sure the prefix is correctly entered.\n", key_fp);
		exit(0);
	}

	sprintf(key_fp, "%s.private_key", key_prefix);

	if (check_file(key_fp))
	{
		fprintf(stderr, "Loading %s\n", key_fp);
		std::ifstream pri_key_file;
		pri_key_file.open(key_fp, std::ifstream::in);
		secret_key.load(context, pri_key_file);
		pri_key_file.close();
	}
	else
	{
		fprintf(stderr, "Private key does not exist @ %s -- Make sure the prefix is correctly entered.\n", key_fp);
		exit(0);
	}

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	bfvHEmpute hempute(encryptor, evaluator, decryptor);

	long nsnp = tag_variant_coord_regs->size();
	long n_test = sample_ids->size();

	FILE* f_decrypted_plaintext_op = open_f(decrypted_geno_op_fp, "w");

	double scale = W0scale;
	std::ifstream f_encrypted_geno_matrix;
	f_encrypted_geno_matrix.open(encrypted_geno_matrix_fp, std::ifstream::in);
	for (int i_tag_snp = 0; i_tag_snp < nsnp; i_tag_snp++)
	{		
		Ciphertext encXData;
		encXData.load(context, f_encrypted_geno_matrix);

		Plaintext plaintmp;
		decryptor.decrypt(encXData, plaintmp);

		uint64_t* plainvec = plaintmp.data();       // convert hex into decimal
		vector<uint64_t> umsg;
		for (long i = 0; i < n_test; ++i) {
			umsg.push_back(plainvec[i]);
		}

		vector<double> dmsg;
		hempute.decode_Vector(dmsg, umsg, plaintext_modulus, scale);
	
		fprintf(stderr, "Reading %d/%d SNP (%d samples)     \r", i_tag_snp, nsnp, dmsg.size());
		fprintf(f_decrypted_plaintext_op, "%s\t%d\t%d\t%s_%d", tag_variant_coord_regs->at(i_tag_snp)->chrom, 
			translate_coord(tag_variant_coord_regs->at(i_tag_snp)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
			translate_coord(tag_variant_coord_regs->at(i_tag_snp)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			tag_variant_coord_regs->at(i_tag_snp)->chrom, 
			tag_variant_coord_regs->at(i_tag_snp)->start);

		for (int i_s = 0; i_s < n_test; i_s++)
		{
			fprintf(f_decrypted_plaintext_op, "\t%.3f", (double)(dmsg[i_s] * W0scale) / Wscale);
		} // i_s loop.
		fprintf(f_decrypted_plaintext_op, "\n");

	} // i_tag_snp loop.
	f_encrypted_geno_matrix.close();

	close_f(f_decrypted_plaintext_op, decrypted_geno_op_fp);
}

// Tags are read from file.
void load_HEmpute_linear_model(char* model_fp, 
								vector<double>* model_intercepts,
								vector<vector<double>*>* model_slopes,
								vector<int>* target_variant_coords,
								vector<vector<int>*>* per_target_tag_coords,
								int& n_vicinity)
{
	fprintf(stderr, "Loading the model from %s\n", model_fp);
	FILE* f_model = open_f(model_fp, "r");
	if (f_model == NULL)
	{
		fprintf(stderr, "Could not open the model file.\n");
		exit(0);
	}

	while (1)
	{
		char* cur_line = getline(f_model);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");

		int n_vicinity_w_slope = toks->size() / 2;
		n_vicinity = (n_vicinity_w_slope - 1) / 2;
		//fprintf(stderr, "N_vicinity: %d\n", n_vicinity);

		vector<double>* per_tag_param = new vector<double>();

		// Read the half and load intercept and slope.
		for (int i_par = 0; i_par < n_vicinity_w_slope; i_par++)
		{
			double cur_model_par = atof(toks->at(i_par)->str());
			if (i_par == 0)
			{
				model_intercepts->push_back(cur_model_par);
			}
			else
			{
				per_tag_param->push_back(cur_model_par);
			}
		} // i_par loop.
		model_slopes->push_back(per_tag_param);

		if (per_tag_param->size() != (2 * n_vicinity))
		{
			fprintf(stderr, "Could not load the parameters: %d\n", per_tag_param->size());
			exit(0);
		}

		vector<int>* cur_target_tag_coords = new vector<int>();
		for (int i_par = n_vicinity_w_slope; 
				i_par < (n_vicinity_w_slope + 2*n_vicinity); 
				i_par++)
		{
			int cur_tag_posn = atoi(toks->at(i_par)->str());
			
			cur_target_tag_coords->push_back(cur_tag_posn);
		} // i_par loop.

		per_target_tag_coords->push_back(cur_target_tag_coords);

		int cur_target_coord = atoi(toks->back()->str());

		target_variant_coords->push_back(cur_target_coord);

		t_string::clean_tokens(toks);
		delete[] cur_line;
	} // param file reading loop.

	fprintf(stderr, "Loaded %d-vicinity models for %d targets.\n", n_vicinity, target_variant_coords->size());
	fclose(f_model);
}


void plaintext_evaluate_LM(char* tag_geno_signals_fp,
	char* sample_ids_list_fp,
	char* model_fp)
{
	fprintf(stderr, "Loading model parameters.\n");
	vector<double>* hempute_model_intercepts = new vector<double>();
	vector<vector<double>*>* hempute_model_slopes = new vector<vector<double>*>();

	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_signal_regions_wrapper(tag_geno_signals_fp, sample_ids_list_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	vector<int>* target_variant_coords = new vector<int>();
	vector<vector<int>*>* per_target_tag_coords = new vector<vector<int>*>();
	int n_vicinity = 0;
	load_HEmpute_linear_model(model_fp, hempute_model_intercepts,
								hempute_model_slopes, target_variant_coords, per_target_tag_coords, n_vicinity);

	fprintf(stderr, "Assigning tags to targets.\n");
	vector<long> tag_model_starting_index;
	int tag_i = 0;
	for (int i_target = 0; i_target < target_variant_coords->size(); i_target++)
	{
		// Move the tag index back.
		while (tag_i > 0 &&
			tag_genotype_signal_regs->at(tag_i)->start >= target_variant_coords->at(i_target))
		{
			tag_i--;
		}

		while (tag_i < tag_genotype_signal_regs->size() &&
			tag_genotype_signal_regs->at(tag_i)->start <= target_variant_coords->at(i_target))
		{
			tag_i++;
		}

		//fprintf(stderr, "Out of loop tag_i: %d (%d) for target %d\n", tag_i, tag_genotype_signal_regs->at(tag_i)->start, target_variant_coords->at(i_target));

		int end_tag_i = tag_i + n_vicinity - 1;
		int start_tag_i = tag_i - n_vicinity;

		if (tag_i < n_vicinity)
		{
			start_tag_i = 0;
		}

		if ((tag_i + n_vicinity - 1) >= tag_genotype_signal_regs->size())
		{
			end_tag_i = tag_genotype_signal_regs->size() - 1;
		}

		// Add a sanity check using the loaded tag coords.
		if (tag_genotype_signal_regs->at(start_tag_i)->start != (per_target_tag_coords->at(i_target)->at(0) + 1))
		{
			fprintf(stderr, "Sanity check failed: target: %d (%d; %d), tag_i: %d (%d); Before-After:[%d-%d]\n",
				i_target, target_variant_coords->at(i_target), per_target_tag_coords->at(i_target)->at(0),
				start_tag_i, tag_genotype_signal_regs->at(start_tag_i)->start,
				tag_genotype_signal_regs->at(start_tag_i - 1)->start, tag_genotype_signal_regs->at(start_tag_i + 1)->start);

			exit(0);
		}

		if (tag_genotype_signal_regs->at(end_tag_i)->start != (per_target_tag_coords->at(i_target)->back() + 1))
		{
			fprintf(stderr, "Sanity check failed: target: %d (%d; %d), tag_i: %d (%d); Before-After:[%d-%d]\n",
				i_target, target_variant_coords->at(i_target), per_target_tag_coords->at(i_target)->back(),
				end_tag_i, tag_genotype_signal_regs->at(end_tag_i)->start,
				tag_genotype_signal_regs->at(end_tag_i - 1)->start, tag_genotype_signal_regs->at(end_tag_i + 1)->start);

			exit(0);
		}

		// Add the starting index to the list.
		tag_model_starting_index.push_back(start_tag_i);
	} // i_target loop.

	FILE* f_plaintext_geno = open_f("plaintext_genotypes.sig", "w");
	for (int i_target = 0; i_target < tag_model_starting_index.size(); i_target++)
	{
		fprintf(f_plaintext_geno, "%s\t%d\t%d\t%s_%d", 
			tag_genotype_signal_regs->at(0)->chrom,
			target_variant_coords->at(i_target),
			target_variant_coords->at(i_target)+1,
			tag_genotype_signal_regs->at(0)->chrom, target_variant_coords->at(i_target));

		for (int i_s = 0; i_s < sample_ids->size(); i_s++)
		{
			double cur_sample_geno_sig = hempute_model_intercepts->at(i_target);
			for (int tag_i = tag_model_starting_index.at(i_target);
				tag_i < tag_model_starting_index.at(i_target) + n_vicinity * 2;
				tag_i++)
			{
				void** cur_tag_info = (void**)(tag_genotype_signal_regs->at(tag_i)->data);
				char* geno_sigs = (char*)(cur_tag_info[0]);

				cur_sample_geno_sig += (double)(geno_sigs[i_s]) * 0.5 * hempute_model_slopes->at(i_target)->at(tag_i - tag_model_starting_index.at(i_target));
			} // tag_i loop.

			fprintf(f_plaintext_geno, "\t%.4f", cur_sample_geno_sig);
		} // i_s loop.

		fprintf(f_plaintext_geno, "\n");
	} // i_target loop.
	fclose(f_plaintext_geno);
}

// Note that this runs on the server.
void secure_evaluate(char* encrypted_tag_geno_fp, 
	char* sample_ids_list_fp,
	char* tag_coords_bed_fp,
	char* public_key_fp,
	char* model_fp,				// Model is trained on the server;; is only specific on the browser.	
	char* enc_target_geno_op_fp,
	char* target_coords_op_bed_fp,
	char* server_update_cmd)
{
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	int n_test_samples = sample_ids->size();

	fprintf(stderr, "Loading model parameters.\n");
	vector<double>* hempute_model_intercepts = new vector<double>();
	vector<vector<double>*>* hempute_model_slopes = new vector<vector<double>*>();

	vector<int>* target_variant_coords = new vector<int>();
	vector<vector<int>*>* per_target_tag_coords = new vector<vector<int>*>();
	int n_vicinity = 0;
	load_HEmpute_linear_model(model_fp, hempute_model_intercepts,
								hempute_model_slopes, target_variant_coords, per_target_tag_coords, n_vicinity);

	vector<t_annot_region*>* tag_genotype_signal_regs = load_BED(tag_coords_bed_fp);

	fprintf(stderr, "Assigning tags to targets.\n");
	vector<long> tag_model_starting_index;
	int tag_i = 0;
	for (int i_target = 0; i_target < target_variant_coords->size(); i_target++)
	{
		// Move the tag index back until tag is to the left of target.
		while (tag_i > 0 &&
				tag_genotype_signal_regs->at(tag_i)->start >= (target_variant_coords->at(i_target)+1))
		{
			tag_i--;
		}

		// Assumes tags and targets do not overlap; there is no point in imputing a variant that already exists.
		while ((tag_i+1) < tag_genotype_signal_regs->size())
		{ 
			if (hempute_OLS_closest_left_tag_checker(tag_genotype_signal_regs, tag_i,
													target_variant_coords->at(i_target)+1))
			{
				break;
			}

			//if (tag_genotype_signal_regs->at(tag_i)->start <= (target_variant_coords->at(i_target) + 1) &&
			//	tag_genotype_signal_regs->at(tag_i + 1)->start > (target_variant_coords->at(i_target) + 1))
			//{
			//	break;
			//}
			
			tag_i++;
		}

		//fprintf(stderr, "Out of loop tag_i: %d (%d) for target %d\n", tag_i, tag_genotype_signal_regs->at(tag_i)->start, target_variant_coords->at(i_target));
		// Tag i is now pointing to the left tag variant on the left side of this target variant.
		int end_tag_i = tag_i + n_vicinity;
		int start_tag_i = tag_i - n_vicinity + 1;

		if (tag_i < n_vicinity)
		{
			start_tag_i = 0;
		}

		if ((tag_i + n_vicinity - 1) >= tag_genotype_signal_regs->size())
		{
			end_tag_i = tag_genotype_signal_regs->size() - 1;
		}

		// Add a sanity check using the loaded tag coords.
		if (tag_genotype_signal_regs->at(start_tag_i)->start != (per_target_tag_coords->at(i_target)->at(0) + 1))
		{
			fprintf(stderr, "Sanity check failed -- start: target: %d (%d; %d), tag_i: %d (%d); Before-After:[%d-%d]\n", 
					i_target, target_variant_coords->at(i_target), per_target_tag_coords->at(i_target)->at(0),
					start_tag_i, tag_genotype_signal_regs->at(start_tag_i)->start,
					tag_genotype_signal_regs->at(start_tag_i - 1)->start, tag_genotype_signal_regs->at(start_tag_i + 1)->start);

			exit(0);
		}

		if (tag_genotype_signal_regs->at(end_tag_i)->start != (per_target_tag_coords->at(i_target)->back() + 1))
		{
			fprintf(stderr, "Sanity check failed -- end: target: %d (%d; %d), tag_i: %d (%d); Before-After:[%d-%d]\n",
				i_target, target_variant_coords->at(i_target), per_target_tag_coords->at(i_target)->back(),
				end_tag_i, tag_genotype_signal_regs->at(end_tag_i)->start,
				tag_genotype_signal_regs->at(end_tag_i - 1)->start, tag_genotype_signal_regs->at(end_tag_i + 1)->start);

			exit(0);
		}

		// Add the starting index to the list.
		tag_model_starting_index.push_back(start_tag_i);
	} // i_target loop.

	long dim = 2 * n_vicinity;

	// Initialize the encryption context.
	EncryptionParameters parms(scheme_type::BFV);
	size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
	int logq = 27;

	// The plain_modulus does not play much of a role;
	// we choose some reasonable value t > results so that it satisfies (results mod t) = results
	size_t plaintext_modulus = (1 << 10);

	double Xscale = 2.0;
	double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
	double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg

	parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { logq }));
	parms.set_plain_modulus(plaintext_modulus);

	auto context = SEALContext::Create(parms);
	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	// Read and load the public key.
	if (check_file(public_key_fp))
	{
		fprintf(stderr, "Loading %s\n", public_key_fp);
		std::ifstream pub_key_file;
		pub_key_file.open(public_key_fp, std::ifstream::in);
		public_key.load(context, pub_key_file);
		pub_key_file.close();
	}
	else
	{
		fprintf(stderr, "Public key does not exist @ %s -- Make sure the prefix is correctly entered.\n", public_key_fp);
		exit(0);
	}

	//sprintf(key_fp, "%s.private_key", key_prefix);

	//if (check_file(key_fp))
	//{
	//	fprintf(stderr, "Loading %s\n", key_fp);
	//	std::ifstream pri_key_file;
	//	pri_key_file.open(key_fp, std::ifstream::in);
	//	secret_key.load(context, pri_key_file);
	//	pri_key_file.close();
	//}
	//else
	//{
	//	fprintf(stderr, "Private key does not exist @ %s -- Make sure the prefix is correctly entered.\n", key_fp);
	//	exit(0);
	//}

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	bfvHEmpute hempute(encryptor, evaluator, decryptor);

	// Start reading the tag genotype matrix.
	std::ifstream f_encrypted_tag_geno_matrix;
	f_encrypted_tag_geno_matrix.open(encrypted_tag_geno_fp, std::ifstream::in);

	/*
	Determine which tags were saved for this model, this will determine what we load here.
	*/
	fprintf(stderr, "Reading %d tag genotypes.\n", tag_genotype_signal_regs->size());
	vector<Ciphertext> enc_tag_geno;
	for (int i_tag_snp = 0; i_tag_snp < tag_genotype_signal_regs->size(); i_tag_snp++)
	{
		fprintf(stderr, "Reading %d/%d SNP         \r", i_tag_snp, tag_genotype_signal_regs->size());
		Ciphertext encXData;
		encXData.load(context, f_encrypted_tag_geno_matrix);

		enc_tag_geno.push_back(encXData);
	} // i_tag_snp loop.

	fprintf(stderr, "\nLoaded encypted tag genotypes.\n");

	f_encrypted_tag_geno_matrix.close();

	// Copy the intercept and slopes.
	int n_targets = hempute_model_intercepts->size();
	fprintf(stderr, "Copying slopes and intercepts for %d targets\n", hempute_model_intercepts->size());
	vector<double> hempute_model_intercept_array;
	hempute_model_intercept_array.resize(n_targets);

	vector<vector<double>> hempute_model_slope_mat;
	hempute_model_slope_mat.resize(n_targets, vector<double>(2*n_vicinity+1));

	for (int i_target = 0; i_target < hempute_model_intercepts->size(); i_target++)
	{
		//hempute_model_intercept_array.push_back(hempute_model_intercepts->at(i_target));
		hempute_model_intercept_array[i_target] = hempute_model_intercepts->at(i_target);

		//vector<double> cur_target_slope_array;
		for (int i_tag = 0; i_tag < hempute_model_slopes->at(i_target)->size(); i_tag++)
		{
			//cur_target_slope_array.push_back(hempute_model_slopes->at(i_target)->at(i_tag));
			hempute_model_slope_mat[i_target][i_tag] = hempute_model_slopes->at(i_target)->at(i_tag);
		} // i_tag loop.

		//hempute_model_slope_mat.push_back(cur_target_slope_array);
	} // i_target loop.
	fprintf(stderr, "Done!\n");

	// Load the encrypted tag data.
	//bfvHEmpute::HEimpute(vector<Ciphertext>& res, vector<bool>& encres_index,
	//	vector<Ciphertext> encXData, dvec model0, dmat model, long n_test,
	//	double W0scale, double Wscale, size_t plaintext_modulus, vector<long> tag_index)
	fprintf(stderr, "Allocating HEimpute object.\n");	
	vector<Ciphertext> encres(n_targets);
	vector<bool> encres_index(n_targets, false);

	vector<Ciphertext> res(n_targets);

	hempute.HEimpute(res, encres_index,
					enc_tag_geno, hempute_model_intercept_array, hempute_model_slope_mat, n_test_samples,
					W0scale, Wscale, plaintext_modulus, tag_model_starting_index);

	// Save the results.
	fprintf(stderr, "Saving %d encrypted imputed target genotypes.\n", res.size());
	ofstream f_enc_target_geno;
	f_enc_target_geno.open(enc_target_geno_op_fp, ios::binary);
	for (int i_target = 0; i_target < res.size(); i_target++)
	{
		res.at(i_target).save(f_enc_target_geno);
	} // i_target loop.
	f_enc_target_geno.close();

	// Write the target variant positions.
	char* chrom_id = tag_genotype_signal_regs->at(0)->chrom;
	FILE* f_target_coords_op = open_f(target_coords_op_bed_fp, "w");
	for (int i_target = 0; i_target < target_variant_coords->size(); i_target++)
	{
		// The weird indexing comes from how the coordinates are written into the model file.
		fprintf(f_target_coords_op, "%s\t%d\t%d\n", 
				chrom_id, 
				target_variant_coords->at(i_target), 
				target_variant_coords->at(i_target) + 1);
	} // i_target loop.
	fclose(f_target_coords_op);
}

// TODO: Get rid of target variant coords for this function, it is not needed here, just save all the tag variants.
// TODO: We don't need vicinity parameter here.
void generic_encrypt_genotypes(char* tag_genotype_signal_matbed_fp,
									char* sample_ids_list_fp, 
									char* key_prefix,
									char* enc_geno_op_fp)
{
	//vector<t_annot_region*>* target_variant_coord_regs = load_BED(target_variant_coords_bed_fp);
	//fprintf(stderr, "Loaded %d target variant coordinates.\n", target_variant_coord_regs->size());

	int n_vicinity = 0;
	long dim = 2 * n_vicinity;

	EncryptionParameters parms(scheme_type::BFV);
	size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
	int logq = 27;

	// The plain_modulus does not play much of a role;
	// we choose some reasonable value t > results so that it satisfies (results mod t) = results
	size_t plaintext_modulus = (1 << 10);

	double Xscale = 2.0;
	double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
	double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg

	parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { logq }));
	parms.set_plain_modulus(plaintext_modulus);

	auto context = SEALContext::Create(parms);
	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	// Write the public/private keys.
	fprintf(stderr, "Loading keys.\n");
	char key_fp[1000];
	sprintf(key_fp, "%s.public_key", key_prefix);

	if (check_file(key_fp))
	{
		fprintf(stderr, "Loading %s\n", key_fp);
		std::ifstream pub_key_file;
		pub_key_file.open(key_fp, std::ifstream::in);
		public_key.load(context, pub_key_file);
		pub_key_file.close();
	}
	else
	{
		fprintf(stderr, "Public key does not exist @ %s -- Make sure the prefix is correctly entered.\n", key_fp);
		exit(0);
	}

	sprintf(key_fp, "%s.private_key", key_prefix);

	if (check_file(key_fp))
	{
		fprintf(stderr, "Loading %s\n", key_fp);
		std::ifstream pri_key_file;
		pri_key_file.open(key_fp, std::ifstream::in);
		secret_key.load(context, pri_key_file);
		pri_key_file.close();
	}
	else
	{
		fprintf(stderr, "Private key does not exist @ %s -- Make sure the prefix is correctly entered.\n", key_fp);
		exit(0);
	}

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	bfvHEmpute hempute(encryptor, evaluator, decryptor);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_signal_regions_wrapper(tag_genotype_signal_matbed_fp, sample_ids_list_fp);
	fprintf(stderr, "Loaded %d genotype signal regions for %d samples.\n", tag_genotype_signal_regs->size(), sample_ids->size());

	vector<vector<int>> tag_geno_data;	

	fprintf(stderr, "Building the tag matrix from %d tag regions.\n", tag_genotype_signal_regs->size());
	for (int i_reg = 0; i_reg < tag_genotype_signal_regs->size(); i_reg++)
	{
		// Generate the genotype signal.
		void** cur_reg_info = (void**)(tag_genotype_signal_regs->at(i_reg)->data);
		char* cur_geno_sig = (char*)(cur_reg_info[0]);

		vector<int> cur_reg_geno_sig;

		for (int i_s = 0; i_s < sample_ids->size(); i_s++)
		{
			cur_reg_geno_sig.push_back(cur_geno_sig[i_s]);
		} // i_s loop.

		tag_geno_data.push_back(cur_reg_geno_sig);
	} // i_reg loop.
	fprintf(stderr, "Built the tag matrix of length %d.\n", tag_geno_data.size());

	fprintf(stderr, "Assigning tags to targets.\n");
	vector<long> tag_model_starting_index;

	// This is the last index that will be used for encryption.
	tag_model_starting_index.push_back(tag_genotype_signal_regs->size());

	//
	fprintf(stderr, "Encrypting..\n");
	vector<Ciphertext> encXData;
	hempute.encrypt_data(encXData, tag_geno_data, tag_model_starting_index, dim);

	// Save the ciphertexts.
	fprintf(stderr, "Saving the encrypted genotypes to %s of length %d\n", enc_geno_op_fp, encXData.size());

	ofstream f_enc_geno;
	f_enc_geno.open(enc_geno_op_fp, ios::binary);
	for (int i_enc = 0; i_enc < encXData.size(); i_enc++)
	{
		encXData.at(i_enc).save(f_enc_geno);
	} // i_enc loop.
	f_enc_geno.close();
}

bool validate_enc_tag_genotype_loading(char* enc_geno_fp,
										char* tag_var_coords_fp,
										char* sample_ids_list_fp)
{
	int n_vicinity = 0;
	long dim = 2 * n_vicinity;

	EncryptionParameters parms(scheme_type::BFV);
	size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
	int logq = 27;

	// The plain_modulus does not play much of a role;
	// we choose some reasonable value t > results so that it satisfies (results mod t) = results
	size_t plaintext_modulus = (1 << 10);

	double Xscale = 2.0;
	double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
	double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg

	parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { logq }));
	parms.set_plain_modulus(plaintext_modulus);

	auto context = SEALContext::Create(parms);
	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	bfvHEmpute hempute(encryptor, evaluator, decryptor);

	std::ifstream f_enc_geno;
	f_enc_geno.open(enc_geno_fp, std::ifstream::in);

	vector<t_annot_region*>* tag_var_coords = load_BED(tag_var_coords_fp);	

	int n_read_vars = 0;
	for (int i_enc = 0; i_enc < tag_var_coords->size(); i_enc++)
	{
		// Read the current encrypted data.
		Ciphertext encXData;
		//if (encXData.load(context, f_enc_geno) == std::invalid_argument ||
		//	encXData.load(context, f_enc_geno) == std::logic_error ||
		//	encXData.load(context, f_enc_geno) == std::runtime_error)
		try
		{
			encXData.load(context, f_enc_geno);
		}
		catch (exception e)
		{
			fprintf(stderr, "Caught exception.\n");
			return(false);
		}
		
		n_read_vars++;
	} // i_enc loop.

	if (f_enc_geno.peek() != istream::traits_type::eof())
	{
		fprintf(stderr, "File not ended after reading all ciphertext.\n");
		return(false);
	}

	f_enc_geno.close();

	fprintf(stderr, "Read %d ciphertexts successfully.\n", n_read_vars);

	return(true);
}

void filter_enc_genotypes(char* enc_geno_fp,
	char* tag_var_coords_fp,
	char* cur_filtering_vars_BED_fp,
	char* sample_ids_list_fp,
	char* filt_enc_geno_fp,
	char* filt_vars_BED_fp)
{
	int n_vicinity = 0;
	long dim = 2 * n_vicinity;

	EncryptionParameters parms(scheme_type::BFV);
	size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
	int logq = 27;

	// The plain_modulus does not play much of a role;
	// we choose some reasonable value t > results so that it satisfies (results mod t) = results
	size_t plaintext_modulus = (1 << 10);

	double Xscale = 2.0;
	double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
	double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg

	parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { logq }));
	parms.set_plain_modulus(plaintext_modulus);

	auto context = SEALContext::Create(parms);
	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	bfvHEmpute hempute(encryptor, evaluator, decryptor);

	// We do not need the sample size for filtering.
	//vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	vector<t_annot_region*>* tag_var_coords = load_BED(tag_var_coords_fp);
	for (int i_reg = 0; i_reg < tag_var_coords->size(); i_reg++)
	{
		tag_var_coords->at(i_reg)->data = NULL;
	} // i_reg loop.

	vector<t_annot_region*>* filt_tag_var_coords = load_BED(cur_filtering_vars_BED_fp);

	// This is necessary to ensure we extract regions at the same place.
	vector<t_annot_region*>* intersects = intersect_annot_regions(tag_var_coords, filt_tag_var_coords, true);
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* tag_var_reg = int_info->src_reg;
		tag_var_reg->data = int_info->dest_reg;
	} // i_int loop.

	vector<t_annot_region*>* check_tag_var_coords = load_BED(tag_var_coords_fp);
	sort(tag_var_coords->begin(), tag_var_coords->end(), sort_regions);
	if (check_tag_var_coords->size() != tag_var_coords->size())
	{
		fprintf(stderr, "The intersected and original tag variants do not match in length: %d/%d\n", 
				check_tag_var_coords->size(), tag_var_coords->size());
		exit(0);
	}

	// This check ensure that encrypted data save protects the order of variants with the original variants.
	for (int i_tag = 0; i_tag < check_tag_var_coords->size(); i_tag++)
	{
		if (check_tag_var_coords->at(i_tag)->start != tag_var_coords->at(i_tag)->start)
		{
			fprintf(stderr, "The intersected and original tag variants do not match in coords: %d/%d\n",
					check_tag_var_coords->at(i_tag)->start, tag_var_coords->at(i_tag)->start);
			exit(0);
		}
	} // i_tag loop.

	std::ifstream f_enc_geno;
	f_enc_geno.open(enc_geno_fp, std::ifstream::in);

	ofstream f_filt_enc_geno;
	f_filt_enc_geno.open(filt_enc_geno_fp, ios::binary);

	vector<t_annot_region*>* dumped_tag_var_regs = new vector<t_annot_region*>();
	for (int i_enc = 0; i_enc < tag_var_coords->size(); i_enc++)
	{
		// Read the current encrypted data.
		Ciphertext encXData;
		encXData.load(context, f_enc_geno);

		if (tag_var_coords->at(i_enc)->data != NULL)
		{
			dumped_tag_var_regs->push_back(tag_var_coords->at(i_enc));
			encXData.save(f_filt_enc_geno);
		}
		else
		{
			// Skip this variant; it was not overlapped.
			fprintf(stderr, "Skipping %s:%d\n", 
				tag_var_coords->at(i_enc)->chrom, 
				tag_var_coords->at(i_enc)->start);
		}
	} // i_enc loop.
	f_enc_geno.close();
	f_filt_enc_geno.close();

	dump_BED(filt_vars_BED_fp, dumped_tag_var_regs);
}

// Divide VCF into chromosomes and save as matbed file.
// TODO: Read chromosome id's from the VCF.
void validate_tag_variants_VCF_write_interim(char* tag_vcf_fp, char* array_platform, int max_n_samples, char* interim_dir)
{
	// Read the VCF file, separate into chromosomes, write chromosome list, write sample id's (anonymized) and save to the intermediate directory.
	vector<char*>* chr_ids = new vector<char*>();

	// Initialize file pointers.
	vector<FILE*>* per_chr_f_ptrs = new vector<FILE*>();
	//for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	//{
	//	per_chr_f_ptrs[i_chr] = NULL;
	//} // i_chr loop.

	//fprintf(stderr, "Processing %d chromosomes.\n", chr_ids->size());

	vector<char*>* anon_tag_sample_ids = NULL;
	vector<char*>* orig_tag_sample_ids = NULL;

	FILE* f_tag_vcf = open_f(tag_vcf_fp, "r");
	int n_vars = 0;
	while (1)
	{
		n_vars++;
		if (n_vars % 10000 == 0)
		{
			fprintf(stderr, "@ %d. tag variant.        \r");
		}

		char* cur_line = getline(f_tag_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		char cur_chrom[100];
		sscanf(cur_line, "%s", cur_chrom);

		if (t_string::compare_strings(cur_chrom, "#CHROM"))
		{
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");
			if (!toks->at(8)->starts_with("FORMAT"))
			{
				fprintf(stderr, "Could not find GT entry in CHROM entry of VCF.\n");
				exit(0);
			}
			else
			{
				anon_tag_sample_ids = new vector<char*>();
				orig_tag_sample_ids = new vector<char*>();

				// Load the tag sample id's.
				for (int i_tok = 9; i_tok < toks->size(); i_tok++)
				{
					char anonymized_sample_id[1000];
					sprintf(anonymized_sample_id, "sample_%d", i_tok-9);
					orig_tag_sample_ids->push_back(t_string::copy_me_str(toks->at(i_tok)->str()));
					anon_tag_sample_ids->push_back(t_string::copy_me_str(anonymized_sample_id));
				} // i_tok loop.

				fprintf(stderr, "Loaded %d sample id's.\n", anon_tag_sample_ids->size());

				if (anon_tag_sample_ids->size() > max_n_samples)
				{
					fprintf(stderr, "HECRYPT can process samples of size less than %d individuals, there are currently %d individuals. Please decrease the sample size.\n", 
							max_n_samples, anon_tag_sample_ids->size());

					exit(0);
				}

				t_string::clean_tokens(toks);
			}
		}
		
		// Skip the comment and sample lines after processing the sample line.
		if (cur_chrom[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if (orig_tag_sample_ids == NULL || 
			anon_tag_sample_ids == NULL)
		{
			fprintf(stderr, "Could not find header line in the VCF file that starts with \"#CHROM\"\n");
			exit(0);
		}

		int i_chr = t_string::get_i_str(chr_ids, cur_chrom);
		if (i_chr == chr_ids->size())
		{
			char cur_vcf_fp[1000];
			sprintf(cur_vcf_fp, "%s/%s.vcf.gz", interim_dir, cur_chrom);
			FILE* f_cur_chr_vcf = open_f(cur_vcf_fp, "w");
			if (f_cur_chr_vcf == NULL)
			{
				fprintf(stderr, "Could not open interim VCF @ %s\n", cur_vcf_fp);
				exit(0);
			}

			chr_ids->push_back(t_string::copy_me_str(cur_chrom));
			per_chr_f_ptrs->push_back(f_cur_chr_vcf);

			i_chr = t_string::get_i_str(chr_ids, cur_chrom);
			if (i_chr == chr_ids->size())
			{
				fprintf(stderr, "Cannot find the chromosome\n");
				exit(0);
			}
		} // chr id existence check.

		// If we have this chromosome, store it.
		fprintf(per_chr_f_ptrs->at(i_chr), "%s\n", cur_line);

		delete [] cur_line;
	} // vcf file reading loop.
	close_f(f_tag_vcf, tag_vcf_fp);

	char sample_ids_fp[1000];
	sprintf(sample_ids_fp, "%s/sample_ids.list", interim_dir);
	FILE* f_sample_ids = open_f(sample_ids_fp, "w");
	if (f_sample_ids == NULL)
	{
		fprintf(stderr, "Could not open sample id's file @ %s\n", interim_dir);
		exit(0);
	}

	for (int i_s = 0; i_s < anon_tag_sample_ids->size(); i_s++)
	{
		fprintf(f_sample_ids, "%s\n", anon_tag_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_sample_ids);

	// Copy the original sample ids.
	sprintf(sample_ids_fp, "original_sample_ids.list");
	f_sample_ids = open_f(sample_ids_fp, "w");
	if (f_sample_ids == NULL)
	{
		fprintf(stderr, "Could not open sample id's file @ %s\n", sample_ids_fp);
		exit(0);
	}

	for (int i_s = 0; i_s < orig_tag_sample_ids->size(); i_s++)
	{
		fprintf(f_sample_ids, "%s\n", orig_tag_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_sample_ids);
	
	// Store the existing chromosome's in the VCF file, not the ones that don't exist.
	char existing_chr_ids_list_fp[1000];
	sprintf(existing_chr_ids_list_fp, "%s/chr_ids.list", interim_dir);
	FILE* f_existing_chr_ids = open_f(existing_chr_ids_list_fp, "w");

	// Binarize the vcf's and delete them.
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		// If this chromosome exists, save it.
		fprintf(f_existing_chr_ids, "%s\n", chr_ids->at(i_chr));

		fprintf(stderr, "Converting %s\n", chr_ids->at(i_chr));
		char cur_vcf_fp[1000];
		sprintf(cur_vcf_fp, "%s/%s.vcf.gz", interim_dir, chr_ids->at(i_chr));
		close_f(per_chr_f_ptrs->at(i_chr), cur_vcf_fp);

		// Save the current tags.
		char tag_reg_bed_fp[1000];
		sprintf(tag_reg_bed_fp, "%s/%s.bed", interim_dir, chr_ids->at(i_chr));
		vector<t_annot_region*>* tag_regs = load_VCF_regions(cur_vcf_fp, false);

		// Get the unique entries: Remove any duplicates.
		vector<t_annot_region*>* unique_tag_regs = get_unique_regions_by_posn(tag_regs, true);
		fprintf(stderr, "Found %d unique regions by position out of %d, excluding repeated variants.\n", unique_tag_regs->size(), tag_regs->size());

		// Sort tag regions before binarizing. This is the only last place that the variants are sorted. The tag regions should not be sorted again,
		// otherwise the imputation will fail since it depends on tracing the tag variants in the vicinity of target variants.
		sort(unique_tag_regs->begin(), unique_tag_regs->end(), sort_regions);
		dump_BED(tag_reg_bed_fp, unique_tag_regs);

		// Binarize this vcf file.
		fprintf(stderr, "Saving variants.\n");
		char cur_bin_geno_fp[1000];
		sprintf(cur_bin_geno_fp, "%s/%s.matbed.gz", interim_dir, chr_ids->at(i_chr));
		extract_genotype_signals_per_VCF(cur_vcf_fp,
			sample_ids_fp,
			tag_reg_bed_fp,
			chr_ids->at(i_chr),
			interim_dir,
			false,
			false,
			false,
			cur_bin_geno_fp);
	} // i_chr loop.

	fclose(f_existing_chr_ids);
}

void generate_keys(char* out_path_prefix)
					//dmat& ypred, 
					//dvec model0, 
					//dmat model, 
					//vector<long> tag_index
{
	chrono::high_resolution_clock::time_point time_start, time_end;
	struct rusage usage;

	//long p = model.size();          // nsnptarget_model = 80,882

	// HE parameters
	EncryptionParameters parms(scheme_type::BFV);
	size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
	int logq = 27;

	// The plain_modulus does not play much of a role;
	// we choose some reasonable value t > results so that it satisfies (results mod t) = results
	size_t plaintext_modulus = (1 << 10);

	double Xscale = 2.0;
	double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
	double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg

	parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { logq }));
	parms.set_plain_modulus(plaintext_modulus);

	cout << "+------------------------------------+" << endl;
	cout << "|           Generating Key           |" << endl;
	cout << "+------------------------------------+" << endl;

	time_start = chrono::high_resolution_clock::now();

	fprintf(stderr, "SEAL Context starting.");
	auto context = SEALContext::Create(parms);
	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	//Encryptor encryptor(context, public_key);
	//Evaluator evaluator(context);
	//Decryptor decryptor(context, secret_key);

	time_end = chrono::high_resolution_clock::now();
	auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	cout << "Scheme generation (milliseconds) : " << time_diff.count() / 1000.0 << endl;
	getrusage(RUSAGE_SELF, &usage);
	cout << "RAM Usage (GB): " << setprecision(4) << (double)usage.ru_maxrss / (DATAParam::memoryscale) << endl;

	// Write the public/private keys.
	char key_fp[1000];
	sprintf(key_fp, "%s.public_key", out_path_prefix);
	std::ofstream pub_key_file;
	pub_key_file.open(key_fp, std::ofstream::out);
	public_key.save(pub_key_file);
	pub_key_file.close();

	sprintf(key_fp, "%s.private_key", out_path_prefix);
	std::ofstream pri_key_file;
	pri_key_file.open(key_fp, std::ofstream::out);
	secret_key.save(pri_key_file);
	pri_key_file.close();
}

