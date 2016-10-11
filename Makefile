###################################
# Makefile written by Brad Nelson #
#           bnelsj@uw.edu         #
###################################

CXX=g++ -lstdc++ -fopenmp -lz -lm

all: bam_chunker_parallel bam_chunker_cascade bam_chunker

bam_chunker_parallel: bin libhts.a
	cd bin && $(CXX) -O3 -I ../htslib/ -L ../htslib -lz ../src/chunker_parallel.cpp ../htslib/libhts.a -lpthread -o bam_chunker_parallel

bam_chunker_cascade: bin libhts.a
	cd bin && $(CXX) -O3 -I ../htslib/ -L ../htslib -lz ../src/chunker_cascade.cpp ../htslib/libhts.a -lpthread -o bam_chunker_cascade

bam_chunker: bin libhts.a
	cd bin && $(CXX) -O3 -I ../htslib/ -L ../htslib -lz ../src/chunker.cpp ../htslib/libhts.a -lpthread -o bam_chunker
bin:
	mkdir bin

libhts.a:
	cd htslib && make
