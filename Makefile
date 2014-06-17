CC := gcc
CPP := g++
LINK := $(CPP)

include make.inc

MODULES := make.bd.inc
# You shouldn't have to go below here

DIRNAME = `dirname $1`
MAKEDEPS = -@gcc -MM -MG $2 $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

# look for include files in each of the modules
MAKEINCLUDES = $(patsubst %, -I%, $1)

# determine the object files
SRCTYPES = c cu cpp

MAKEOBJS = $(foreach srctype, $(SRCTYPES), $(patsubst %.$(srctype), obj/%.o, $(wildcard $(patsubst %, %/*.$(srctype), ($1)))))

# include the C include dependencies
MAKEDEPSINC = $(patsubst obj/%.o, .dep/%.d, $1)

.DEFAULT_GOAL := all

INCLUDES :=
DEPS :=

include $(MODULES)

INCLUDES := $(strip $(patsubst %, %\n, $(INCLUDES)))
INCLUDES := $(shell printf " $(INCLUDES)\n" | sort | uniq)

CFLAGS += $(INCLUDEFLAGS) $(INCLUDES)
CPPFLAGS += $(INCLUDEFLAGS) $(INCLUDES)
NVCCFLAGS += $(INCLUDEFLAGS) $(INCLUDES)

.PHONY : all

# calculate include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CPPFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

.dep/%.d : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

obj/%.o : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CC) $(CFLAGS) -c -o $@ $<

ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif
