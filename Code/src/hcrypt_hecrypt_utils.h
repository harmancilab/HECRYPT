#ifndef __HECRYPT_UTILS__
#define __HECRYPT_UTILS__

void validate_tag_variants_VCF_write_interim(char* tag_vcf_fp, char* array_platform, int max_n_samples, char* interim_dir);

void generic_encrypt_genotypes(char* tag_genotype_signal_matbed_fp, 
								char* sample_ids_list_fp,
								char* key_prefix,
								char* enc_geno_op_fp);

void generic_decrypt_genotypes(char* encrypted_geno_matrix_fp,
								char* geno_coords_fp,
								char* sample_ids_list_fp,
								char* key_prefix,
								const char* op_fp);

void get_R2_per_imputed_genotypes(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp, char* stats_op_fp);

bool validate_enc_tag_genotype_loading(char* enc_geno_fp,
	char* tag_var_coords_fp,
	char* sample_ids_list_fp);

void filter_enc_genotypes(char* enc_geno_fp,
	char* tag_var_coords_fp,
	char* cur_filtering_vars_BED_fp,
	char* sample_ids_list_fp,
	char* filt_enc_geno_fp,
	char* filt_vars_BED_fp);

void plaintext_evaluate_LM(char* tag_geno_signals_fp,
	char* sample_ids_list_fp,
	char* model_fp);

void compare_LMSE_vs_secure_eval(char* lmse_sig_fp, char* secure_val_sig_fp, char* sample_ids_list_fp);

void secure_evaluate(char* encrypted_tag_geno_fp,
	char* sample_ids_list_fp,
	char* tag_coords_bed_fp,
	char* public_key_fp,
	char* model_fp,				// Model is trained on the server;; is only specific on the browser.	
	char* enc_target_geno_op_fp,
	char* target_coords_op_bed_fp,
	char* server_update_cmd);

void generate_keys(char* out_path_prefix);

#endif // __HECRYPT_UTILS__