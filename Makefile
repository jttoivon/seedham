package = seedham
version = 0.9.2
tarname = $(package)
distdir = $(tarname)-$(version)
prefix=/usr/local
datarootdir=$(prefix)/share
docdir=$(datarootdir)/doc/$(package)

ifdef BOOSTROOT
	BOOSTINC=-I $(BOOSTROOT)/include
	BOOSTLIB=-L $(BOOSTROOT)/lib
	PROGOPT=$(BOOSTROOT)/lib/libboost_program_options.a
else
	BOOSTINC=
	BOOSTLIB=
	PROGOPT=-lboost_program_options -lboost_system -lboost_filesystem
endif

CXXFLAGS= -std=gnu++11 -DPACKAGE_VERSION=\"$(version)\" -Wall -Wno-sign-compare -g $(BOOSTINC)

NOOPENMP?=0
ifeq ($(NOOPENMP),0)
	CXXFLAGS += -fopenmp #-Wno-unknown-pragmas
else
	CXXFLAGS += -Wno-unknown-pragmas
endif

DEBUG?=0
OBJDIR=.
ifeq ($(DEBUG),1)
#    CFLAGS =-g3 -gdwarf2 -DDEBUG
	CXXFLAGS +=-O0
	OBJDIR=debug
else
ifeq ($(DEBUG),2)
	CXXFLAGS+=-O0 -pg
	OBJDIR=debug
else
	CXXFLAGS+=-O3
endif
endif

#OPTIMIZE=-O3

# -I /usr/include/x86_64-linux-gnu/c++/4.7/

PROGRAMS=multinomial

CXX=g++
#CXX=$(HOME)/usr/bin/g++

LDFLAGS=$(BOOSTLIB) $(PROGOPT) #-L$(HOME)/usr/lib #-B/usr/bin/ld.gold

all: CPM03 ahocorasick $(addprefix $(OBJDIR)/, $(PROGRAMS))

$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: CPM03
CPM03:
	@$(MAKE) -C CPM03 all

.PHONY: ahocorasick
ahocorasick/libahocorasick.a:
	@$(MAKE) -C ahocorasick

test: test/test_suffix_array_wrapper


install: all
	install -d $(prefix)/bin
	install -d $(docdir)
	install -m 0755 multinomial $(prefix)/bin
	install -m 0755 seedham $(prefix)/bin
	install -m 0755 seedham+ $(prefix)/bin
	install -m 0644 README $(docdir)

dist: $(distdir).tar.gz

$(distdir).tar.gz: FORCE $(distdir)
	tar chf - $(distdir) | gzip -9 -c > $(distdir).tar.gz
	rm -rf $(distdir)

FORCE:
	-rm $(distdir).tar.gz &> /dev/null
	-rm -rf $(distdir) &> /dev/null

distcheck: $(distdir).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(distdir) all check clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz\
          ready for distribution."

check: all
#	./multinomial --multimer 9 data/TFAP2A-head-1000.seq > /dev/null
	./seedham 8 data/TFAP2A-head-1000.seq > /dev/null
	./seedham+ 8 data/TFAP2A-head-1000.seq > /dev/null
	@echo "*** ALL TESTS PASSED ***"

$(distdir):
	mkdir -p $(distdir)
	mkdir -p $(distdir)/CPM03
	mkdir -p $(distdir)/data	
	mkdir -p $(distdir)/ahocorasick	
	mkdir -p $(distdir)/ahocorasick/doc	
	cp seedham $(distdir)
	cp seedham+ $(distdir)
	cp README $(distdir)
	cp COPYING $(distdir)
	cp Makefile $(distdir)
	cp common.cpp $(distdir)
	cp probabilities.cpp $(distdir)
	cp parameters.cpp $(distdir)
	cp matrix_tools.cpp $(distdir)
	cp my_assert.cpp $(distdir)
	cp combinatorics.cpp $(distdir)
	cp multinomial_helper.cpp $(distdir)
	cp bndm.cpp $(distdir)
	cp orientation.cpp $(distdir)
	cp data.cpp $(distdir)
	cp iupac.cpp $(distdir)
	cp suffix_array_wrapper.cpp $(distdir)
	cp iupac.hpp $(distdir)
	cp matrix.hpp $(distdir)
	cp vectors.hpp $(distdir)
	cp timing.hpp $(distdir)
	cp matrix_tools.hpp $(distdir)
	cp common.hpp $(distdir)
	cp data.hpp $(distdir)
	cp cob.hpp $(distdir)
	cp orientation.hpp $(distdir)
	cp probabilities.hpp $(distdir)
	cp parameters.hpp $(distdir)
	cp combinatorics.hpp $(distdir)
	cp my_assert.hpp $(distdir)
	cp multinomial_helper.hpp $(distdir)
	cp suffix_array_wrapper.hpp $(distdir)
	cp bndm.hpp $(distdir)
	cp type.hpp $(distdir)
	cp CPM03/checker.hpp $(distdir)/CPM03
	cp CPM03/COPYING $(distdir)/CPM03
	cp CPM03/difference_cover.cpp $(distdir)/CPM03
	cp CPM03/difference_cover.hpp $(distdir)/CPM03
	cp CPM03/doubling.hpp $(distdir)/CPM03
	cp CPM03/Makefile $(distdir)/CPM03
	cp CPM03/partition.hpp $(distdir)/CPM03
	cp CPM03/README $(distdir)/CPM03
	cp CPM03/stringsort.hpp $(distdir)/CPM03
	cp CPM03/suffixsort.hpp $(distdir)/CPM03
	cp CPM03/test-suffixsort.cpp $(distdir)/CPM03
	cp CPM03/timing.hpp $(distdir)/CPM03
	cp data/TFAP2A-head-1000.seq $(distdir)/data
	cp aho_corasick_wrapper.hpp $(distdir)
	cp lambda.cpp		    $(distdir)
	cp lambda.hpp		    $(distdir)
	cp kmp.cpp		    $(distdir)
	cp kmp.hpp		    $(distdir)
	cp dependence.cpp	    $(distdir)
	cp dependence.hpp	    $(distdir)
	cp packed_string.cpp	    $(distdir)
	cp packed_string.hpp	    $(distdir)
	cp kmer_tools.cpp	    $(distdir)
	cp kmer_tools.hpp	    $(distdir)
	cp unordered_map.hpp	    $(distdir)
	cp multinomial.cpp	    $(distdir)
	cp seed_basic.cpp	    $(distdir)
	cp seed_basic.hpp	    $(distdir)
	cp huddinge.cpp	    $(distdir)
	cp huddinge.hpp	    $(distdir)
	cp ahocorasick/ac_types.h     $(distdir)/ahocorasick
	cp ahocorasick/aho_corasick.c $(distdir)/ahocorasick
	cp ahocorasick/aho_corasick.h $(distdir)/ahocorasick
	cp ahocorasick/config.h	      $(distdir)/ahocorasick
	cp ahocorasick/Makefile	      $(distdir)/ahocorasick
	cp ahocorasick/node.c	      $(distdir)/ahocorasick
	cp ahocorasick/node.h	      $(distdir)/ahocorasick
	cp ahocorasick/README	      $(distdir)/ahocorasick
	cp ahocorasick/doc/how_to_use	      $(distdir)/ahocorasick/doc



# programs


MULTINOMIAL_OBJS=matrix_tools.o common.o parameters.o kmp.o data.o\
	dependence.o bndm.o probabilities.o iupac.o\
	my_assert.o combinatorics.o lambda.o kmer_tools.o multinomial_helper.o packed_string.o  \
	suffix_array_wrapper.o orientation.o huddinge.o multinomial.o seed_basic.o
$(OBJDIR)/multinomial: $(addprefix $(OBJDIR)/, $(MULTINOMIAL_OBJS)) ahocorasick/libahocorasick.a CPM03/difference_cover.o
	$(CXX) $(CXXFLAGS) $+ -o $@ $(LDFLAGS)



# test programs

TEST_SUFFIX_ARRAY_WRAPPER_OBJS=CPM03/difference_cover.o suffix_array_wrapper.o iupac.o common.o test/test_suffix_array_wrapper.o
test/test_suffix_array_wrapper: $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS) -o test/test_suffix_array_wrapper $(LDFLAGS)

.PHONY: dist FORCE distcheck


########################
#
# compilation units
#
########################

main_units=multinomial.o

other_units=permutation_test.o dependence.o matrix_tools.o kmp.o bndm.o cob.o vectors.o stats.o cooccurrence.o cooccurrence_utils.o common.o esko.o sum_method.o alignment.o parameters.o bernoulli.o dynamic_programming_common.o probabilities.o combinatorics.o lambda.o my_assert.o kmer_tools.o nucleosome.o orientation.o multinomial_helper.o huddinge.o packed_string.o data.o iupac.o suffix_array_wrapper.o


# Automatically create dependency files (*.d)
%.d: %.cpp
	@set -e; rm -f $@; \
         $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
         sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
         rm -f $@.$$$$

-include $(other_units:.o=.d) $(main_units:.o=.d)


$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@





.PHONY: clean
clean:
	rm -f *.o $(PROGRAMS)
	@$(MAKE) -C CPM03 clean
	@$(MAKE) -C ahocorasick clean

.PHONY: showprograms
showprograms:
	@echo $(PROGRAMS) | tr ' ' '\n'
