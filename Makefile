###################################
# Makefile written by Brad Nelson #
#           bnelsj@uw.edu         #
###################################

CXX=g++ -lstdc++

bam_chunker: bin libhts.a
	cd bin && $(CXX) -O3 -I ../htslib/ -L ../htslib -lz ../src/chunker.cpp ../htslib/libhts.a -lpthread -o bam_chunker
bin:
	mkdir bin
libhts.a:
	cd htslib && make
