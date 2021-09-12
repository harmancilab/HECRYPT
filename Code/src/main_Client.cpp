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
				--array: This option is currently unused.\n\
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
		-clean_encrypted_dir [Directory name]\n\
			Cleans plaintext genotype information written in preprocessing.\n\
	Decryption:\n\
		-decrypt_target_genotypes --results [Encrypted results directory] --key_prefix [Key prefix] --out [Output directory]\n\
			Decrypts and saves the final imputation results.\n\
			Arguments\n\
				--results: Path for the encrypted directory of results that is downloaded from imputation server.\n\
				--key_prefix: Path for the private key. This file has extension \".private_key\".\n\
	Genotype Accuracy:\n\
		-get_accuracy_statistics --imputed_dir [Imputed genotypes directory] --known_dir [Known genotypes directory]\n", argv[0]);

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
	else
	{
		fprintf(stderr, "Unknown option: %s\n", argv[1]);
		exit(0);
	}
}