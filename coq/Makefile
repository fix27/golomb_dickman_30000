# Требуется coq_makefile из состава Coq.
COQ_MAKEFILE ?= coq_makefile
COQ_PROJECT  ?= _CoqProject

all:
	$(COQ_MAKEFILE) -f $(COQ_PROJECT) -o Makefile.coq
	$(MAKE) -f Makefile.coq

clean:
	-$(MAKE) -f Makefile.coq clean
	-$(RM) Makefile.coq Makefile.coq.conf

.PHONY: all clean
