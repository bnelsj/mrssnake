/*

This program was created at:  Fri Feb 26 14:47:00 2016
This program was created by:  Brad Nelson


Contact: bnelsj@gmail.com

Organization: University of Washington

The MIT License (MIT)

Copyright (c) <2016> <Brad Nelson>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"

#define HTS_FMT_BAI 1

struct interval
{
    int64_t start,end;
	bool operator < (const interval& other) const
	{
		return (start < other.start);
	}
};

struct options{
   std::string file;
   int part;
   int nparts;
   int chunk_size;
}globalOpts;

static const char *optString = "b:p:n:c";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    globalOpts.file = "NA";
	globalOpts.chunk_size = 36;
    opt = getopt(argc, argv, optString);
    while(opt != -1){
	switch(opt){
		case 'h':
		{
		 std::cerr << "Useage" << std::endl;
		}
		case 'b':
		{
		 globalOpts.file = optarg;
		 break;
		}
		case 'p':
		{
		 globalOpts.part = atoi(optarg);
		 break;
		}
		case 'n':
		{
		 globalOpts.nparts = atoi(optarg);
		 break;
		}
		case 'c':
		{
		 globalOpts.chunk_size = atoi(optarg);
		 break;
		}
		case '?':
		{
		 break;
		}
	}
 
  opt = getopt( argc, argv, optString ); 
   }
return 1;
}

inline void initKstring(kstring_t * k){
  k->m = 0;
  k->l = 0;
  k->s = 0;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : nchunks in bam, partition number, total number of partitions

 Function does   : calculates start and end chunk to parse

 Function returns: half-open chunk interval for given partition

*/

std::vector<struct interval> get_chunk_range(std::vector<struct interval> chunks, int part, int nparts) {
	struct interval chunk_range;
	int range_size;
	int nchunks = chunks.size();

	range_size = chunks.size() / nparts;
	chunk_range.start = range_size * part;
	part == nparts-1 ? chunk_range.end = nchunks : chunk_range.end = range_size * part + range_size;
	
	//chunk_range.end = range_size * part + range_size;
	
	std::vector<struct interval> chunks_to_read(&chunks[chunk_range.start], &chunks[chunk_range.end]);

	return chunks_to_read;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : Vector of virtual offset intervals, bamfile

 Function does   : Iterate over vector, seek to start of interval and read to end

 Function returns: 0 if successful, 1 otherwise

*/
int get_reads(std::vector<struct interval> &chunks, std::string bamfile)
{
	BGZF* bam;
	
	kstring_t tmp, seq;
	initKstring(&tmp);
	initKstring(&seq);
	int seek;
	if(bgzf_is_bgzf == 0) {
		std::cerr << "Error: file is not in BGZF format." << std::endl;
		exit(1);
	}
	bam = bgzf_open(bamfile.c_str(), "r");
	// Get BGZF file handle

    for(int i = 0; i < chunks.size(); i++){
		int state = 1;
		seek = bgzf_seek(bam, chunks[i].start, 0);
		if(seek != 0) {
			std::cerr << "Error: could not seek to position " << chunks[i].start << "." << std::endl;
			exit(1);
		}
		std::cout << "Going to position " << chunks[i].start << " " << chunks[i].end << std::endl;

		while(state > 0 and state < chunks[i].end) {
		//	for(int j =0; j < 9; j++) {
		//		state = bgzf_getline(bam, '\t', &tmp);
		//	}
			state = bgzf_getline(bam, '\n', &seq);
		//	state = bgzf_getline(bam, '\n', &tmp);

			std::cout << seq.s << std::endl;
		}
		
		//Seek to offset
		//Read until end
			//Split reads
    }
	bgzf_close(bam);

    return 0;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : options

 Function does   : get number of blocks in file, calculate start and end block from part and nparts

 Function returns: start and end positions

*/

std::vector<struct interval> read_index(std::string fn, int part, int nparts) {
	const char * file_name = fn.c_str();
	std::string index_name = fn + ".bai";
	long chunk_counter=0;
	FILE *ptr_myfile;
	std::vector <struct interval> chunks, chunks_to_read;

	ptr_myfile = fopen(index_name.c_str(),"rb");

	if(!ptr_myfile) {
		std::cerr << "Unable to open file." << std::endl;
		exit(1);
	}

	char magic[4];
	fread(&magic, sizeof(char), 4, ptr_myfile);	
	printf("magic:%s\n",magic);

    int32_t n_ref;
    fread(&n_ref,sizeof(int32_t),1,ptr_myfile);
    printf("n_ref:%d\n",n_ref);

    int32_t n_bin;
    fread(&n_bin,sizeof(int32_t),1,ptr_myfile);
    printf("n_bin:%d\n",n_bin);

	int32_t i;
    for (i = 0; i < n_bin; ++i) {
        uint32_t bin; 
        fread(&bin,sizeof(uint32_t),1,ptr_myfile);

        int32_t n_chunk;
        fread(&n_chunk,sizeof(int32_t),1,ptr_myfile);
		chunk_counter += n_chunk;
        int32_t j;
        for (j = 0; j < n_chunk; ++j) {
            int64_t chunk_beg, chunk_end;
			struct interval tmp;
            fread(&chunk_beg,sizeof(int64_t),1,ptr_myfile);
            fread(&chunk_end,sizeof(int64_t),1,ptr_myfile);

			tmp.start = chunk_beg;
			tmp.end = chunk_end;
			chunks.push_back(tmp);		
        }

    }
	printf("Total chunks: %d\n", chunk_counter);
	chunks_to_read = get_chunk_range(chunks, part, nparts);	
	fclose(ptr_myfile);
	printf("Vector size: %d\n", chunks.size());
	printf("Chunks to read: %d\n", chunks_to_read.size());


	return chunks_to_read;
	
}



//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
/*
std::string chunk_read(std::string &read, int chunk_size)
{
	std::stringstream stream;
	int n_to_do = read->length / chunk_size;

    for(i = 0; i < n_to_do; i++){
        stream << ">0" << std::endl;
        stream << read->seq[i*chunk_size : i*chunk_size + chunk_size] << std::endl;
    }

    return stream;
}
*/
//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{

    globalOpts.chunk_size = 36;
	globalOpts.part = -1;
	globalOpts.nparts = -1;
	std::vector<struct interval> chunks_to_read;
    int parse = parseOpts(argc, argv);

    if (globalOpts.part == -1) {
        std::cerr << "Error: no partition specified" << std::endl;
        exit(1);
    }

     if (globalOpts.nparts == -1) {
        std::cerr << "Error: number of partitions not specified" << std::endl;
        exit(1);
    }

	if (!(globalOpts.part < globalOpts.nparts)) {
		std::cerr << "Error: partition number must be less than total number of partitions (Partition numbers are 0-based)" << std::endl;
		exit(1);
	}

    if (globalOpts.file.empty()) {
        std::cerr << "Error: no bam file provided" << std::endl;
		exit(1);
    }

	samFile *in = sam_open(globalOpts.file.c_str(), "r");
	if(in == NULL) {
    	std::cerr << "Unable to open BAM/SAM file." << std::endl;
		exit(1);
    }
	hts_idx_t *idx = sam_index_load(in, globalOpts.file.c_str());

	chunks_to_read = read_index(globalOpts.file, globalOpts.part, globalOpts.nparts);
	int result = get_reads(chunks_to_read, globalOpts.file.c_str());	

    //    idx = get_positions(globalOpts);
//    std::cout << idx << std::endl;   
// Get index from bamfile
// Open index
// Get total chunks
// Calculate start_pos and end_pos
// Get offset of start_pos and end_pos
// Seek to start_pos
   // Chunk reads
// Read until end_pos

    return 0;
}
