#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include "hcrypt_hecrypt_utils.h"
#include "hcrypt_ansi_cli.h"
#include "hcrypt_config.h"
#include "hcrypt_ansi_string.h"
#include "hcrypt_annot_region_tools.h"
#include "hcrypt_variation_tools.h"
#include "hcrypt_utils.h"

#include <vector>
using namespace std;

#include <algorithm>

// cmake for seal:: cmake.org/files/v3.12/cmake-3.12.3.tar.gz

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [option] [arguments]\n\
Options:\n\
	VCF Processing:\n\
		-preprocess_tags_genotypes --VCF [Tag Variant VCF Path] --array [Array platform] --interim [Directory Path for intermediate files]\n\
			Reads and validates a VCF file that contains the input genotypes. This option \n\
			also writes the intermediate files that will be used for encryption only. These\n\
			files can be deleted after encryption.\n\
			Arguments:\n\
				--VCF: The VCF (version >4.0) formatted file that \n\
					contains the tag genotype variants. (<1000 samples)\n\
				--array: The array platform. Currently only Illumina Duo v3 is supported.\n\
				--interim: The directory where intermediate files will be saved for encryption.\n\
					This directory can be deleted after encryption.\n\
	Encryption:\n\
		-generate_key_pair --key_prefix [Prefix for key pair]\n\
			Generates the public-private key pair to encrypt/decrypt genotype matrix. The keys are\n\
			stored in an output file. After the keys are generated, you can copy the private key to a\n\
			confidential location. Encryption requires only public key.\n\
			Arguments:\n\
				--key_prefix: Prefix for key pair that will be appended to the output.\n\
		-encrypt_tag_genotypes --key_prefix [Key file prefix] --interim_dir [Intermediate file]\n\
			Arguments:\n\
				--key_prefix: Prefix for key pair that will be appended to the output.\n\
				--interim_dir: Output encrypted genotype directory\n\
		-validate_encrypted_genotypes --enc_geno_dir [Encrypted genotypes directory]\n\
		-clean_encrypted_dir [Directory name]\n\
			Cleans plaintext genotype information written in preprocessing.\n\
	Evaluation:\n\
		-filter_encrypted_genotypes --filtering_vars [Variants BED file path] --enc_geno_dir [Interim directory] --op_dir [Output directory]\n\
			Filters the encrypted genotypes per variant list in a BED file. Does not decrypt the data.\n\
			Arguments:\n\
				--filtering_vars: BED file that contains the variants to select\n\
				--enc_geno_dir: Encrypted genotypes data directory\n\
				--op_dir: Directory that will be used to save filtered encrypted genotypes\n\
		-secure_evaluate --test_samples [Test samples list file] --enc_geno [Encrypted genotypes file path] --tag_coords [Tag coordinates (Client variants)]\n\
			Secure imputation using the encrypted tag genotypes as input. The test sample id's is an \n\
			anonymized list of samples. Also, the tag coordinates of the genotype matrix must be supplied.\n\
	Decryption:\n\
		-decrypt_target_genotypes --results [Encrypted results directory] --key_prefix [Key prefix] --out [Output directory]\n\
			Decrypts and saves the final imputation results.\n\
			Arguments\n\
				--results: Path for the encrypted directory of results that is downloaded from imputation server.\n\
				--key_prefix: Path for the private key. This file has extension \".private_key\".\n\
	Genotype Accuracy:\n\
		-get_accuracy_statistics --imputed_dir [Imputed genotypes directory] --known_dir [Known genotypes directory]\n\
	Diagnostics:\n\
		-compare_LMSE_vs_secure_eval\n\
		-compare_extracted_genotypes\n\
		-plaintext_impute\n", argv[0]);

		exit(0);
	}

	if (strcmp(argv[1], "-clean_encrypted_dir") == 0)
	{
		if (argc < 3)
		{
			fprintf(stderr, "USAGE: %s %s [Encrypted genotypes directory]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		char* dir_name = argv[2];

		fprintf(stderr, "Continue with cleaning directory \"%s\"? [Y/N]", dir_name);
		if (getc(stdin) != 'y' && getc(stdin) != 'Y')
		{
			fprintf(stderr, "Aborted cleanup.\n");
			exit(0);
		}
		
		fprintf(stderr, "Deleting VCFs..\n");

		char rm_cmd[1000];
		sprintf(rm_cmd, "rm -f %s/*.matbed.gz %s/*.vcf.gz %s/*.vcf", dir_name, dir_name, dir_name);
		system(rm_cmd);

		// Clean interim directory:
		fprintf(stderr, "\
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\
Done! Yon can compress the interim directory:\n\
>> tar -cvjf %s.tar.bz2 %s\n\
After compression is finished, you can start\n\
uploading the encrypted data to it to\n\
secureomics.org/Web.\n\
Please make sure your private key and folder id\n\
are stored in an easy-to-find and secure place\n\
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n", dir_name, dir_name);
	} // -clean_encrypted_dir
	else if (strcmp(argv[1], "-plaintext_impute") == 0)
	{
		if (argc < 5)
		{
			fprintf(stderr, "USAGE: %s %s [Tag genotypes signal regions BED] [Sample id's file path] [Model file path]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		char* tag_geno_signals_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* model_fp = argv[4];

		plaintext_evaluate_LM(tag_geno_signals_fp,
			sample_ids_list_fp,
			model_fp);
	} // -plaintext_impute option.
	else if (strcmp(argv[1], "-get_accuracy_statistics") == 0)
	{
		if (argc < 6)
		{
			fprintf(stderr, "USAGE: %s %s --imputed_dir [Imputed genotypes directory] --known_dir [Known genotypes directory]\n", argv[0], argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;

		char* imputed_dir = t_string::copy_me_str(cli->get_value_by_option("--imputed_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the output directory with \"--imputed_dir\"\n");
			exit(0);
		}

		char* known_dir = t_string::copy_me_str(cli->get_value_by_option("--known_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the model parameters directory with \"--known_dir\"\n");
			exit(0);
		}

		char imp_samples_list_fp[1000];
		sprintf(imp_samples_list_fp, "%s/sample_ids.list", imputed_dir);
		if (!check_file(imp_samples_list_fp))
		{
			fprintf(stderr, "Could not find imputed sample id's @ %s\n", imp_samples_list_fp);
			exit(0);
		}

		char known_samples_list_fp[1000];
		sprintf(known_samples_list_fp, "%s/sample_ids.list", known_dir);
		if (!check_file(known_samples_list_fp))
		{
			fprintf(stderr, "Could not find known sample id's @ %s\n", known_samples_list_fp);
			exit(0);
		}

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", imputed_dir);

		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

		// Process the chromosomes in order.		
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			char cur_chr_imp_geno_fp[1000];
			sprintf(cur_chr_imp_geno_fp, "%s/%s.imp", imputed_dir, chr_ids->at(i_chr));

			char cur_chr_known_geno_fp[1000];
			sprintf(cur_chr_known_geno_fp, "%s/%s.matbed", known_dir, chr_ids->at(i_chr));

			char stats_fp[1000];
			sprintf(stats_fp, "R2_stats_%s.txt", chr_ids->at(i_chr));
			get_R2_per_imputed_genotypes(cur_chr_imp_geno_fp, imp_samples_list_fp,
										cur_chr_known_geno_fp, known_samples_list_fp, stats_fp);
		} // i_chr loop.
	} // -get_accuracy_statistics option.
	else if (strcmp(argv[1], "-compare_LMSE_vs_secure_eval") == 0)
	{
		if (argc < 5)
		{
			fprintf(stderr, "USAGE: %s %s [LMSE signal file path] [Secure imputation signal file path] [Sample id's file path]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		char* lmse_sig_fp = argv[2];
		char* secure_val_sig_fp = argv[3];
		char* sample_ids_list_fp = argv[4];

		compare_LMSE_vs_secure_eval(lmse_sig_fp, secure_val_sig_fp, sample_ids_list_fp);
	} // -compare_LMSE_vs_secure_eval option.
	else if (strcmp(argv[1], "-compare_extracted_genotypes") == 0)
	{
		if (argc < 4)
		{
			fprintf(stderr, "USAGE: %s %s [Original signal matrix] [Extracted signal matrix]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		char* orig_fp = argv[2];
		char* ext_fp = argv[3];

		vector<char*>* orig_lines = buffer_file(orig_fp);
		vector<char*>* ext_lines = buffer_file(ext_fp);

		fprintf(stderr, "Loaded %d, %d regions.\n", orig_lines->size(), ext_lines->size());

		sort(ext_lines->begin(), ext_lines->end(), t_string::sort_strings);

		vector<char*>* matching_lines = new vector<char*>();
		for (int i_orig = 0; i_orig < orig_lines->size(); i_orig++)
		{
			int i_ext = t_string::fast_search_string(orig_lines->at(i_orig), ext_lines, 0, ext_lines->size());

			while (i_ext < ext_lines->size() &&
					i_ext > 0 &&
					(t_string::sort_strings(orig_lines->at(i_orig), ext_lines->at(i_ext)) ||
					t_string::compare_strings(ext_lines->at(i_ext), orig_lines->at(i_orig))))
			{
				i_ext--;
			} // i_ext loop.

			while (i_ext < ext_lines->size() &&
					(t_string::sort_strings(ext_lines->at(i_ext), orig_lines->at(i_orig)) ||
					t_string::compare_strings(ext_lines->at(i_ext), orig_lines->at(i_orig))))
			{
				if (t_string::compare_strings(ext_lines->at(i_ext), orig_lines->at(i_orig)))
				{
					matching_lines->push_back(ext_lines->at(i_ext));
					break;
				}

				i_ext++;
			} // i_ext loop.
		} // i_orig loop.

		fprintf(stderr, "Found %d matching lines.\n", matching_lines->size());
	} // -compare_LMSE_vs_secure_eval option.
	else if (strcmp(argv[1], "-secure_evaluate") == 0)
	{
		if (argc < 8)
		{
			fprintf(stderr, "USAGE: %s %s --model_dir [Model directory] --enc_tag_dir [Encrypted tag genotypes dir] --op_dir [Output directory]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;
		char* encrypted_tag_geno_dir = t_string::copy_me_str(cli->get_value_by_option("--enc_tag_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the encrypted tag directory with \"--enc_tag_dir\"\n");
			exit(0);
		}

		char* op_dir = t_string::copy_me_str(cli->get_value_by_option("--op_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the output directory with \"--op_dir\"\n");
			exit(0);
		}

		char* model_dir = t_string::copy_me_str(cli->get_value_by_option("--model_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the model parameters directory with \"--model_dir\"\n");
			exit(0);
		}

		char test_samples_list_fp[1000];
		sprintf(test_samples_list_fp, "%s/sample_ids.list", encrypted_tag_geno_dir);
		if (!check_file(test_samples_list_fp))
		{
			fprintf(stderr, "Could not find sample id's under encrypted tag genotypes directory.\n");
			exit(0);
		}

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", encrypted_tag_geno_dir);

		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		char public_key_fp[1000];
		sprintf(public_key_fp, "%s/pubkey", encrypted_tag_geno_dir);
		if (!check_file(public_key_fp))
		{
			fprintf(stderr, "Could not find public key @ %s\n", public_key_fp);
			exit(0);
		}

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

		char* server_update_cmd = NULL;

		// Process the chromosomes in order.		
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Processing %s\n",  chr_ids->at(i_chr));

			char model_fp[1000];
			sprintf(model_fp, "%s/%s.model", model_dir, chr_ids->at(i_chr));

			char encrypted_tag_geno_fp[1000];
			sprintf(encrypted_tag_geno_fp, "%s/%s.enc", encrypted_tag_geno_dir, chr_ids->at(i_chr));

			char tag_coords_bed_fp[1000];
			sprintf(tag_coords_bed_fp, "%s/%s.bed", encrypted_tag_geno_dir, chr_ids->at(i_chr));

			char target_coords_op_bed_fp[1000];
			sprintf(target_coords_op_bed_fp, "%s/%s.bed", op_dir, chr_ids->at(i_chr));

			char enc_target_geno_op_fp[1000];
			sprintf(enc_target_geno_op_fp, "%s/%s.enc", op_dir, chr_ids->at(i_chr));
	
			secure_evaluate(encrypted_tag_geno_fp,
							test_samples_list_fp,
							tag_coords_bed_fp,
							public_key_fp,
							model_fp,
							enc_target_geno_op_fp,
							target_coords_op_bed_fp,
							server_update_cmd);
		} // i_chr loop.

		// Copy the chromosome id's and sample id's.
		char imp_chr_ids_fp[1000];
		sprintf(imp_chr_ids_fp, "%s/chr_ids.list", op_dir);
		copy_file(chr_ids_list_fp, imp_chr_ids_fp);

		char imp_sample_ids_list_fp[1000];
		sprintf(imp_sample_ids_list_fp, "%s/sample_ids.list", op_dir);
		copy_file(test_samples_list_fp, imp_sample_ids_list_fp);
	} // -secure_evaluate option.
	else if (strcmp(argv[1], "-decrypt_target_genotypes") == 0)
	{
		if (argc < 8)
		{
			fprintf(stderr, "USAGE: %s %s --key_prefix [Key file prefix] --enc_geno [Encrypted genotype directory] --op_dir [Decrypted genotypes directory]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;
		char* key_prefix = t_string::copy_me_str(cli->get_value_by_option("--key_prefix", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the key prefix with \"--key_prefix\"\n");
			exit(0);
		}

		char* encrypted_geno_dir = t_string::copy_me_str(cli->get_value_by_option("--enc_geno", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the encrypted tag genotypes path with \"--enc_geno\"\n");
			exit(0);
		}

		char* output_dir = t_string::copy_me_str(cli->get_value_by_option("--op_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the output directory with \"--op_dir\"\n");
			exit(0);
		}

		char sample_ids_list_fp[1000];
		sprintf(sample_ids_list_fp, "%s/sample_ids.list", encrypted_geno_dir);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", encrypted_geno_dir);
		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

		// Write each chromosome separately.
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			char encrypted_geno_matrix_fp[1000];
			sprintf(encrypted_geno_matrix_fp, "%s/%s.enc", encrypted_geno_dir, chr_ids->at(i_chr));

			char tag_target_geno_coords_fp[1000];
			sprintf(tag_target_geno_coords_fp, "%s/%s.bed", encrypted_geno_dir, chr_ids->at(i_chr));

			char output_geno_matrix_fp[1000];
			sprintf(output_geno_matrix_fp, "%s/%s.imp", output_dir, chr_ids->at(i_chr));

			generic_decrypt_genotypes(encrypted_geno_matrix_fp,
				tag_target_geno_coords_fp,
				sample_ids_list_fp,
				key_prefix,
				output_geno_matrix_fp);
		} // i_chr loop.

		// Copy the chromosome id's to decrypted directory.
		fprintf(stderr, "Copying chromosome id's and sample id's.\n");
		char dec_chr_ids_fp[1000];
		sprintf(dec_chr_ids_fp, "%s/chr_ids.list", output_dir);
		copy_file(chr_ids_list_fp, dec_chr_ids_fp);

		char orig_sample_fp[1000];
		sprintf(orig_sample_fp, "original_sample_ids.list");
		char dec_sample_fp[1000];
		sprintf(dec_sample_fp, "%s/sample_ids.list", output_dir);
		copy_file(orig_sample_fp, dec_sample_fp);
	} // -decrypt_tag_variants
	else if (strcmp(argv[1], "-encrypt_tag_genotypes") == 0)
	{
		if (argc < 6)
		{
			fprintf(stderr, "USAGE: %s %s --key_prefix [Key file prefix] --interim_dir [Interim directory]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;
		char* key_prefix = t_string::copy_me_str(cli->get_value_by_option("--key_prefix", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the key prefix with \"--key_prefix\"\n");
			exit(0);
		}

		char* interim_dir = t_string::copy_me_str(cli->get_value_by_option("--interim_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the intermediate directory with \"--interim_dir\"\n");
			exit(0);
		}

		long n_vicinity = 0;

		char sample_ids_list_fp[1000];
		sprintf(sample_ids_list_fp, "%s/sample_ids.list", interim_dir);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", interim_dir);		

		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());
		
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Encrypting %s..\n", chr_ids->at(i_chr));
			char enc_geno_op_fp[1000];
			sprintf(enc_geno_op_fp, "%s/%s.enc", interim_dir, chr_ids->at(i_chr));

			char tag_genotype_signal_matbed_fp[1000];
			sprintf(tag_genotype_signal_matbed_fp, "%s/%s.matbed.gz", interim_dir, chr_ids->at(i_chr));

			generic_encrypt_genotypes(tag_genotype_signal_matbed_fp, 
										sample_ids_list_fp,
										key_prefix,
										enc_geno_op_fp);
		} // i_chr loop.

		fprintf(stderr, "Copying public key.\n");
		char pubkey_fp[1000];
		sprintf(pubkey_fp, "%s.public_key", key_prefix);
		char dest_pubkey[1000];
		sprintf(dest_pubkey, "%s/pubkey", interim_dir);
		copy_file(pubkey_fp, dest_pubkey);
	} // -encrypt_tag_variants option.	
	else if (strcmp(argv[1], "-validate_encrypted_genotypes") == 0)
	{
		if (argc < 4)
		{
			fprintf(stderr, "USAGE: %s %s -validate_encrypted_genotypes --enc_geno_dir [Encrypted genotypes directory]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;
		char* enc_geno_dir = t_string::copy_me_str(cli->get_value_by_option("--enc_geno_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the intermediate directory with \"--enc_geno_dir\"\n");
			exit(0);
		}

		char sample_ids_list_fp[1000];
		sprintf(sample_ids_list_fp, "%s/sample_ids.list", enc_geno_dir);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", enc_geno_dir);

		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Testing encrypted genotype data on %s..\n", chr_ids->at(i_chr));

			// Load the existing encrypted genotype data.
			char enc_geno_fp[1000];
			sprintf(enc_geno_fp, "%s/%s.enc", enc_geno_dir, chr_ids->at(i_chr));

			char tag_var_coords_fp[1000];
			sprintf(tag_var_coords_fp, "%s/%s.bed", enc_geno_dir, chr_ids->at(i_chr));

			if (validate_enc_tag_genotype_loading(enc_geno_fp,
													tag_var_coords_fp,
													sample_ids_list_fp) == false)
			{
				fprintf(stderr, "Loading failed.\n");
				FILE* f_stat = open_f("valid.stats", "w");
				fprintf(f_stat, "FAIL");
				fclose(f_stat);
				return 1;
			}
		} // i_chr loop.

		FILE* f_stat = open_f("valid.stats", "w");
		fprintf(f_stat, "SUCCESS");
		fclose(f_stat);
	} // validate_encrypted_genotypes
	else if (strcmp(argv[1], "-filter_encrypted_genotypes") == 0)
	{
		if (argc < 6)
		{
			fprintf(stderr, "USAGE: %s %s --filtering_vars [Variants BED file path] --enc_geno_dir [Interim directory] --op_dir [Output directory]\n",
				argv[0],
				argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;
		char* filter_vars_BED_fp = t_string::copy_me_str(cli->get_value_by_option("--filtering_vars", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the variant BED file path with \"--filtering_vars\"\n");
			exit(0);
		}

		char* enc_geno_dir = t_string::copy_me_str(cli->get_value_by_option("--enc_geno_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the intermediate directory with \"--enc_geno_dir\"\n");
			exit(0);
		}

		char* op_dir = t_string::copy_me_str(cli->get_value_by_option("--op_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the intermediate directory with \"--op_dir\"\n");
			exit(0);
		}

		long n_vicinity = 0;

		char sample_ids_list_fp[1000];
		sprintf(sample_ids_list_fp, "%s/sample_ids.list", enc_geno_dir);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", enc_geno_dir);

		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<t_annot_region*>* filtering_var_regs = load_BED(filter_vars_BED_fp);
		fprintf(stderr, "Loaded %d variants from %s\n", filtering_var_regs->size(), filter_vars_BED_fp);

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

		// Copy sample id's.
		char geno_sample_ids_list_fp[1000];
		sprintf(geno_sample_ids_list_fp, "%s/sample_ids.list", enc_geno_dir);
		char filt_geno_sample_ids_list_fp[1000];
		sprintf(filt_geno_sample_ids_list_fp, "%s/sample_ids.list", op_dir);
		copy_file(geno_sample_ids_list_fp, filt_geno_sample_ids_list_fp);

		// Start processing each chromosome.
		char filtered_chr_ids_fp[1000];
		sprintf(filtered_chr_ids_fp, "%s/chr_ids.list", op_dir);
		FILE* f_filtered_chr_ids = open_f(filtered_chr_ids_fp, "w");
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Filtering %s..\n", chr_ids->at(i_chr));
			vector<t_annot_region*>* cur_chr_filter_var_regs = get_regions_per_chromosome(filtering_var_regs, chr_ids->at(i_chr));
			fprintf(stderr, "%d regions.\n", cur_chr_filter_var_regs->size());

			// Load the existing encrypted genotype data.
			char enc_geno_fp[1000];
			sprintf(enc_geno_fp, "%s/%s.enc", enc_geno_dir, chr_ids->at(i_chr));

			char filt_enc_geno_fp[1000];
			sprintf(filt_enc_geno_fp, "%s/%s.enc", op_dir, chr_ids->at(i_chr));

			// Check if there there are overlaps on this chromosome, this intersection does not do anything, just gets the number of intersects.
			char tag_var_coords_fp[1000];
			sprintf(tag_var_coords_fp, "%s/%s.bed", enc_geno_dir, chr_ids->at(i_chr));
			vector<t_annot_region*>* cur_chr_tag_var_regs = load_BED(tag_var_coords_fp);

			vector<t_annot_region*>* intersects = intersect_annot_regions(cur_chr_filter_var_regs, cur_chr_tag_var_regs, false);
			fprintf(stderr, "Filtering encrypted %d tag genotypes on %s.\n", intersects->size(), chr_ids->at(i_chr));

			// If there are intersects, process the overlaps on this chromosome.
			if (intersects->size() > 0)
			{
				// Write the chromosome id.
				fprintf(f_filtered_chr_ids, "%s\n", chr_ids->at(i_chr));

				char cur_filtering_vars_BED_fp[1000];
				sprintf(cur_filtering_vars_BED_fp, "%s.bed", chr_ids->at(i_chr));
				dump_BED(cur_filtering_vars_BED_fp, cur_chr_filter_var_regs);

				char filt_vars_BED_fp[1000];
				sprintf(filt_vars_BED_fp, "%s/%s.bed", op_dir, chr_ids->at(i_chr));

				// Filter the chromosomes.
				filter_enc_genotypes(enc_geno_fp,
					tag_var_coords_fp,
					cur_filtering_vars_BED_fp,
					sample_ids_list_fp,
					filt_enc_geno_fp,
					filt_vars_BED_fp);
			}
		} // i_chr loop.

		fprintf(stderr, "Copying public key.\n");
		char pubkey_fp[1000];
		sprintf(pubkey_fp, "%s/pubkey", enc_geno_dir);
		if (!check_file(pubkey_fp))
		{
			fprintf(stderr, "Could not find the public key @ %s\n", pubkey_fp);
			exit(0);
		}

		char dest_pubkey[1000];
		sprintf(dest_pubkey, "%s/pubkey", op_dir);
		copy_file(pubkey_fp, dest_pubkey);
	} // filter_encrypted_genotypes
	else if (strcmp(argv[1], "-generate_key_pair") == 0)
	{
		if (argc < 3)
		{
			fprintf(stderr, "USAGE: %s %s --key_prefix [VCF Path]\n", 
					argv[0],
					argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;
		char* key_prefix = t_string::copy_me_str(cli->get_value_by_option("--key_prefix", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the key prefix with \"--key_prefix\"\n");
			exit(0);
		}

		generate_keys(key_prefix);
	} // -generate_key_pair option.
	else if (strcmp(argv[1], "-preprocess_tags_genotypes") == 0)
	{
		if (argc < 8)
		{
			fprintf(stderr, "USAGE: %s %s --VCF [VCF Path] --array [Array platform] --interim [Directory Path for intermediate files]\n", 
					argv[0], 
					argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");

		bool res = false;
		char* vcf_fp = t_string::copy_me_str(cli->get_value_by_option("--VCF", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the VCF file path with \"--VCF\"\n");
			exit(0);
		}

		char* array_platform = t_string::copy_me_str(cli->get_value_by_option("--array", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify array platform with \"--array\"\n");
			exit(0);
		}

		char* interim_dir = t_string::copy_me_str(cli->get_value_by_option("--interim", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify intermediate directory with \"--interim\"\n");
			exit(0);
		}

		if (!check_file(vcf_fp))
		{
			fprintf(stderr, "Could not find VCF file @ %s\n", vcf_fp);
			exit(0);
		}

		// Check interim directory.
		char sample_ids_fp[1000];
		sprintf(sample_ids_fp, "%s/sample_ids.list", interim_dir);
		FILE* f_sample_ids = open_f(sample_ids_fp, "w");
		if (f_sample_ids == NULL)
		{
			fprintf(stderr, "Could not open interim directory directory @ %s\n", interim_dir);
			exit(0);
		}
		else
		{
			fclose(f_sample_ids);
		}

		validate_tag_variants_VCF_write_interim(vcf_fp, array_platform, 1000, interim_dir);


	} // -validate_VCF option.
}