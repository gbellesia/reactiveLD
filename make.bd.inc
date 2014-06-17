BD_RUN_MODULES := src/ include

BD_RUN_OBJS := $(call MAKEOBJS, $(BD_RUN_MODULES))
BD_DEPS := $(call MAKEDEPSINC, $(BD_RUN_OBJS))

bd : bd_run

# link the program
bd_run : $(BD_RUN_OBJS)
	$(LINK) $(LFLAGS) -o $@ $^ $(LIBS)

bd_clean:
	-@rm -f $(BD_RUN_OBJS) $(BD_DEPS) bd_run

INCLUDES += $(call MAKEINCLUDES, $(BD_RUN_MODULES))
DEPS += $(BD_DEPS)
clean : bd_clean
.PHONY : bd_clean bd
all : bd
