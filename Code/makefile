all: HECRYPT_client.bin

CC = g++
gsl_flags = -lz -lpthread -lgsl -lgslcblas 
exec_name = bin/HECRYPT_client.bin
LIB_DIR = src
UTMSR_DIR = ../external/UTMSR/src
SEAL_INCLUDE_DIR = ../external/UTMSR/external/seal/native/src
SEAL_LIB_DIR = ../external/UTMSR/external/seal/native/lib/libseal-3.4.a

comp_flags = -c -O3 -Wall -std=c++17 -I${SEAL_INCLUDE_DIR}

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

%.o: %.c
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
src/main_Client.o \
src/hcrypt_hecrypt_utils.o \
${UTMSR_DIR}/param.o \
${UTMSR_DIR}/utils_data.o \
${UTMSR_DIR}/thread.o \
${UTMSR_DIR}/utils.o \
${UTMSR_DIR}/HEmpute_bfv.o \
${UTMSR_DIR}/HEmpute_ckks.o \
${LIB_DIR}/hcrypt_genome_sequence_tools.o \
${LIB_DIR}/hcrypt_histogram.o \
${LIB_DIR}/hcrypt_xlog_math.o \
${LIB_DIR}/hcrypt_mapped_read_tools.o \
${LIB_DIR}/hcrypt_config.o \
${LIB_DIR}/hcrypt_ansi_cli.o \
${LIB_DIR}/hcrypt_seed_manager.o \
${LIB_DIR}/hcrypt_ansi_thread.o \
${LIB_DIR}/hcrypt_annot_region_tools.o \
${LIB_DIR}/hcrypt_signal_track_tools.o \
${LIB_DIR}/hcrypt_ansi_string.o \
${LIB_DIR}/hcrypt_gff_utils.o \
${LIB_DIR}/hcrypt_utils.o \
${LIB_DIR}/hcrypt_nomenclature.o \
${LIB_DIR}/hcrypt_variation_tools.o \
${LIB_DIR}/hcrypt_rng.o \
${LIB_DIR}/hcrypt_nucleotide.o

HECRYPT_client.bin: ${objs}	
	${CC} -O3 ${gsl_flags} -I ${SEAL_INCLUDE_DIR} -o ${exec_name} ${objs} ${SEAL_LIB_DIR}

clean:
	rm -f ${objs} ${exec_name} 

