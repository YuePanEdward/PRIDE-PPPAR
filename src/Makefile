
SUBDIR = lib \
		 orbit \
		 lsq \
		 arsig \
		 tedit \
		 redig \
		 get_ctrl \
		 utils

all : $(SUBDIR)

$(SUBDIR) :
	make -j4 -C $@

install :
	@mkdir -p ../bin/
	cp arsig/arsig  lsq/lsq tedit/tedit redig/redig orbit/sp3orb \
		orbit/mergesp3 get_ctrl/get_ctrl utils/xyz2enu ../bin/

.PHONY : all $(SUBDIR) install

clean :
	-rm -rf $(addsuffix /.obj, $(SUBDIR))

