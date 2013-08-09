#
# mex/src/misc/Makefile
#

ifndef src_top
	src_top := $(HOME)/l/src
endif

ifndef Makedef
	Makedef := $(src_top)/make/Makedef.$(Arch)
endif

include $(Makedef)

incs=	-I../def \
	$(MatlabInc) -DMATLAB_HOST -DMmex

CFLAGS= $(DEBUG) $(BCFLAGS) $(incs)

doth=	../def/*.h $(MatlabH)

goal1=	$(bin_mex)/delaysum1_mex.mex$(MatlabSuf)

#
# trick: to avoid use of the IO libraries, we disable CountAlloc here
# so that the standard calloc is used.
#
com= $(MatlabCom) -UCountAlloc -v

all:	$(goal1)

g1c = ../def/mexarg.c delaysum1_mex.c
$(goal1):	$(doth) $(g1c) $(libs)
		$(com) -output $@ $(incs) -DNeed_delaysum1_mex_gateway \
			$(g1c) $(libs)
		chmod 755 $@

reset:
	rm -i $(goal1)
