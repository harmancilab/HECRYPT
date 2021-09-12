# HECRYPT: 

This repository contains the documentation for HECRYPT -- genotype encryption/decryption tool for starting imputation with secure imputation server. HECRYPT is a command line tool that runs on Linux systems.

This repository also contains the documentation and example data (under "data/") including VCF data, encrypted data, and the encryption/decryption keys. 

## Build ##

## Standalone Executable ##
Alternatively, you can download HECRYPT from [here](https://secureomics.org/Main/Web/HECRYPT.bin). HECRYPT requires gzip, tar executables; and gsl, zlib libraries to be installed.

After downloading HECRYPT, you need to set it as an executable:
```
chmod 755 HECRYPT.bin
```

### Docker Usage
This executable may fail on certain installations of Linux. In this case, an alternative is to use the docker image. 

You can run the docker container as following:
```
docker pull secureomics/hecrypt:v5
docker run -v $PWD:/host -i -t secureomics/hecrypt:v5 /bin/bash

# gzip/bzip2, gsl and tar are needed to compress/decompress data:
yum -y install gzip bzip2 tar gsl
```

When the container starts, you should be able to see that HECRYPT.bin is located under root directory /. You can copy your genotype data into a running container. For this, open a new command line (docker run command is running as above) and run following:
```
# Following returns the list of running containers, choose the container id.
docker container ls

# Use the container id in the command below:
docker cp data/tag_data.vcf.gz [container id]:/
```
You can customize where the data is copied. After copying the genotype data, you can continue using the container to process and encrypt it.

**Please make sure to make a copy of the private key while working within docker container.**

### Example Dataset

We processed 3000 tag SNVs on chromosomes 18, 19, 20 and generated an example VCF file is included under "data/" directory. 

In this directory there are the intermediate encrypted data (intermediate.tar.bz2) that can be uploaded to server directly. Also, the public/private keys are included for decrypting the results.

## Preprocessing gVCF File that contains the Tag Genotypes ##
HECRYPT accepts a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) formatted file that contains the variants and sample genotypes (starting from the 10th column of the file). This VCF file contains the genotypes of the tag variants that are used in predicting missing target variants in the population panel.

It is necessary to ensure the VCF file is formatted correctly. HECRYPT requires a header line that starts with "#CHROM", otherwise it will not process the VCF file.

To process a VCF file, we run following commands on command line:
```
rm -f -r intermediate
mkdir intermediate
./HECRYPT.bin -preprocess_tags_genotypes --VCF data/tag_data.vcf.gz --array Illumina --interim intermediate
```
After running the preprocessing command, HECRYPT reads and separates the VCF file into chromosomes, and converts them to a format that can be quickly loaded. These intermediate results are written under the directory that is specific by "--interim" option. 

In case HECRYPT complains that there are more than 1000 individuals, you can filter the first 1000 individuals:
```
gzip -cd tag_data.vcf.gz | cut -f1-1009 | gzip > tag_data_smaller.vcf.gz
```
Use the previous command to process the tag genotypes in 'tag_data_smaller.vcf.gz'.

The repeated entries are removed from the VCF file. This is necessary to exclude the redundant tag genotypes or multi-allelic variants, which are not reliably used in imputation models, yet.

**Currently, HECRYPT can process at most 1000 samples in the VCF file. If there are more than 1000 samples, HECRYPT will write an error message and exit.**

## Key Generation:
After VCF is preprocessed, we generate the public/private keys:
```
./HECRYPT.bin -generate_key_pair --key_prefix my_key
```
This command generates and saves two files: 
1. "my_key.public_key": Contains the public key, which is needed for encryption of the tag genotype data.
2. "my_key.private_key": Contains the private key, which is needed for decryption of the imputed variant genotypes.

You can change the prefix of these keys using the "--key_prefix" option of "-generate_key_pair" option.

## Encryption of Genotypes:
After keys are generated, we encrypt the genotypes:
```
./HECRYPT.bin -encrypt_tag_genotypes --key_prefix my_key --interim_dir intermediate
```
Here, "--key_prefix" and "--interim_dir" must match the key prefix from key generation and the name of the directory that we stored the preprocessed data.

The encryption writes the encrypted genotypes to the same directory as the intermediate directory.

## Data Cleaning
The intermediate directory needs to be cleaned before submission to remove the unencrypted 
genotype information that is written while preprocessing. For this, we run:
```
./HECRYPT.bin -clean_encrypted_dir intermediate
```
This command deletes the intermediate files that contain plaintext information. The remaining files include:
1. Sample id's (anonymized as "sample_[#]"), 
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

After this, you can navigate to https://secureomics.org/OpenImpute and start uploading the file named "intermediate.tar.bz2". Note that the name of the file is not relevant since it is never sent to the server. However, server keeps track of a simple hash of this file for resuming interrupted file uploads. This is done to ensure that a different file is not uploaded after a failed upload attempt is being resumed.

## Downloading data 
After imputation is submitted, the server performs imputation. When the results are available for download, a button is displayed where the results can be downloaded from. 
The downloaded file has the name "RESULTS.tar.bz2" by default.

## Decryption of imputed target variant genotypes
After downloading the imputed genotypes, we run following command:
```
tar -xvjf RESULTS.tar.bz2
mkdir decrypted_genotypes
./HECRYPT.bin -decrypt_target_genotypes --key_prefix my_key --enc_geno analysis/RESULTS --op_dir decrypted_genotypes
```
This command decrypts and write the results under "decrypted_genotypes/" for each chromosome separately, where each of them are named "1.imp, 2.imp, ..., 22.imp".

Each file is a text file that contains the soft genotypes that are predicted by imputation server. These values can be rounded upto nearest integer in [0,1,2] to generate final imputed values.