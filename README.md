# HECRYPT: 

This repository contains the documentation for HECRYPT -- genotype encryption/decryption tool for secure imputation. HECRYPT is a command line tool that runs on Linux systems.

## Build ##
You can download HECRYPT from [here](https://secureomics.org/Web/./HECRYPT.bin). HECRYPT requires zlib libraries to be installed, including gzip so that compressed files can be used as input.

After downloading HECRPT, you need to set it as an executable:
```
chmod 755 HECRPT.bin
```

An example VCF file is included under "data/" directory.

## Preprocessing gVCF File that contains the Tag Genotypes ##
HECRYPT accepts a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) formatted file that contains the variants and sample genotypes (starting from the 10th column of the file). This VCF file contains the genotypes of the tag variants that are used in predicting missing target variants in the population panel.

It is necessary to ensure the VCF file is formatted correctly. HECRYPT requires a header line that starts with "#CHROM", otherwise it will not process the VCF file.

To process a VCF file, we run following commands on command line:
```
rm -f -r intermediate
mkdir intermediate
./HECRYPT.bin -preprocess_tags_genotypes --VCF data/example.vcf.gz --array Illumina --interim intermediate
```
After running the preprocessing command, HECRYPT reads and separates the VCF file into chromosomes, and converts them to a format that can be quickly loaded. These intermediate results are written under the directory that is specific by "--interim" option. 
Also, the repeated entries are removed from the 

## Key Generation:
After VCF is preprocessed, we generate the public/private keys:
```
HECRYPT -generate_key_pair --key_prefix my_key
```
This command gnenerates and saves two files: 
1. "my_key.public_key": Contains the public key, which is needed for encryption of the tag genotype data.
2. "my_key.public_key": Contains the private key, which is needed for decryption of the imputed variant genotyes.

You can change the prefix of these keys using the "--key_prefix" option of "-generate_key_pair" option.

## Encryption of Genotypes:
After keys are generated, we encrypt the genotypes:
```
HECRYPT -encrypt_tag_genotypes --key_prefix my_key --interim_dir intermediate
```
Here, "--key_prefix" and "--interim_dir" must match the key prefix from key generation and the name of the directory that we stored the preprocessed data.

The encryption writes the encrypted genotypes to the same directory as the intermediate directory.

## Data Cleaning
The intermediate directory needs to be cleaned before submission to remove the unencrypted 
genotype information that is written while preprocessing. For this, we run:
```
HECRYPT -clean_encrypted_dir intermediate
```
This command deletes the intermediate files that contain plaintext information. The remaining files include:
1. Sample id's (anonymized), 
2. Tag coordinates, 
3. Chromosome identifiers list, 
4. Encrypted genotypes, named in the format [chr_id].enc, 
5. Public key (Cannot be used to decrypt the data).

## Submission to the Server
Final step is compression of the files:
```
tar -cvjf intermediate.tar.bz2 intermediate
```
After starting upload, we recommended to copy the keys to a safe place. We also recommend saving your online folder identifier with the private key.

After this, you can navigate to https://secureomics.org/Web and start uploading the file named "intermediate.tar.bz2". Note that the name of the directory is not relevant since it is never sent to the server. However, server keeps track of a simple hash of this file for resuming interrupted file uploads. This is done to ensure that a different file is not uploaded after a failed upload attempt is being resumed.

## Downloading data 
After imputation is submitted, the server performs imputation. After the results are available for download, a button is displayed where the results can be downloaded from. 
The downloaded file has the named "RESULTS.tar.bz2" by default.

## Decryption of imputed target variant genotypes
After downloading the imputed genotypes, we run following command:
```
tar -xvjf RESULTS.tar.bz2
mkdir decrypted_genotypes
HECRYPT -decrypt_target_genotypes --key_prefix my_key --enc_geno analysis/RESULTS --op_dir decrypted_genotypes
```
This command decrypts and write the results under "decrypted_genotypes/" for each chromosome separately, where each of them are named "1.imp, 2.imp, ..., 22.imp".

Each file is a text file that contains the soft genotypes that are predicted by imputation server. These values can be rounded upto nearest integer in [0,1,2] to generate final imputed values.