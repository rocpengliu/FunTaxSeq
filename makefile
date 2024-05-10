CC = gcc
CXX = g++
CFLAGS = -O3 -DNDEBUG -fno-omit-frame-pointer
CXXFLAGS = -O3 -pthread -std=c++11 -DNDEBUG -g -fno-inline -fno-omit-frame-pointer
LDLIBS = -ldl -lpthread -lz 
#CFLAGS = -O3 -DNDEBUG -fsanitize=address -fno-omit-frame-pointer
#CXXFLAGS = -O3 -pthread -std=c++11 -DNDEBUG -g -fno-inline -fno-omit-frame-pointer -fsanitize=address
#LDLIBS = -ldl -lpthread -lz -fsanitize=address
INCLUDES	= -I./src/include -I./src/include/ncbi-blast+ -I./src/zlib

BLASTOBJS = src/include/ncbi-blast+/algo/blast/core/pattern.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_posit.o \
			 src/include/ncbi-blast+/algo/blast/composition_adjustment/matrix_frequency_data.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_dynarray.o \
			 src/include/ncbi-blast+/algo/blast/core/matrix_freq_ratios.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_encoding.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_stat.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_filter.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_util.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_message.o \
			 src/include/ncbi-blast+/algo/blast/core/ncbi_erf.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_options.o \
			 src/include/ncbi-blast+/algo/blast/core/ncbi_math.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_program.o \
			 src/include/ncbi-blast+/algo/blast/core/ncbi_std.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_psi_priv.o \
			 src/include/ncbi-blast+/util/tables/raw_scoremat.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_query_info.o \
			 src/include/ncbi-blast+/algo/blast/core/blast_seg.o 

BWTOBJS = src/bwt/bwt.o src/bwt/compactfmi.o src/bwt/sequence.o src/bwt/suffixArray.o

ifeq ($(uname -s), "Darwin")
LD_LIBS_STATIC = -Wl,-all_load -lpthread -lz -Wl,-noall_load
else
LD_LIBS_STATIC = -Wl,--whole-archive -lpthread -lz -Wl,--no-whole-archive
endif

all: makefile  src/fts src/ftd src/bwt/mkbwt 
	mkdir -p bin
	cp src/fts bin/
	cp src/ftd bin/
	cp src/bwt/mkbwt bin/mkbwt
	cp src/bwt/mkfmi bin/mkfmi

src/bwt/mkbwt:
	$(MAKE) -C src/bwt/ $(MAKECMDGOALS)
	
src/fts: makefile src/bwt/mkbwt src/fts.o src/homosearcher.o src/transsearcher.o src/dnasearcher.o src/fragment.o src/bwtfmiDB.o src/adaptertrimmer.o src/basecorrector.o \
	src/duplicate.o src/evaluator.o src/fastareader.o src/fastqreader.o src/filter.o src/filterresult.o src/htmlreporter.o \
	src/jsonreporter.o  src/nucleotidetree.o src/options.o src/overlapanalysis.o src/peprocessor.o \
	src/polyx.o src/processor.o src/read.o src/seprocessor.o src/sequence.o src/stats.o src/threadconfig.o src/umiprocessor.o \
	src/unittest.o src/writer.o src/writerthread.o $(BLASTOBJS)
	
	$(CXX) $(LDFLAGS) -o src/fts src/fts.o src/homosearcher.o src/transsearcher.o src/dnasearcher.o src/fragment.o src/bwtfmiDB.o src/adaptertrimmer.o src/basecorrector.o \
	src/duplicate.o src/evaluator.o src/fastareader.o src/fastqreader.o src/filter.o src/filterresult.o src/htmlreporter.o \
	src/jsonreporter.o  src/nucleotidetree.o src/options.o src/overlapanalysis.o src/peprocessor.o \
	src/polyx.o src/processor.o src/read.o src/seprocessor.o src/sequence.o src/stats.o src/threadconfig.o src/umiprocessor.o \
	src/unittest.o src/writer.o src/writerthread.o $(BWTOBJS) $(BLASTOBJS) $(LDLIBS)

src/ftd: makefile src/ftd.o src/phylotree.o src/funtaxoptions.o src/funtaxdecoder.o 
	$(CXX) $(LDFLAGS) -o src/ftd src/ftd.o src/phylotree.o src/funtaxoptions.o src/funtaxdecoder.o $(LDLIBS)

#%.o : %.c makefile
%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
#%.o : %.cpp makefile
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f -v src/bwt/mkbwt src/bwt/mkfmi fts bin/* testdata/All* testdata/*.html testdata/*_mapped* testdata/D*.txt testdata/*.json testdata/*.txt.gz
	find src/ -name "*.o" -delete
	$(MAKE) -C src/bwt/ clean

static: LDFLAGS = -static
static: LDLIBS = $(LD_LIBS_STATIC)
static: all

debug: CXXFLAGS = -O3 -pthread -std=c++11 -g -Wall -Wpedantic -Wextra -Wconversion -fno-omit-frame-pointer
debug: CFLAGS = -g -O3 -Wall -Wno-uninitialized
#debug: CXXFLAGS = -O3 -pthread -std=c++11 -g -Wall -Wpedantic -Wextra -Wconversion -fno-omit-frame-pointer -fsanitize=address
#debug: CFLAGS = -g -O3 -Wall -Wno-uninitialized -fsanitize=address
debug: all

.PHONY: clean debug static