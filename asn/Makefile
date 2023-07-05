CC=mpicc
CFLAGS=-std=c11 # 2017-08: the latest standard
CFLAGS+=-W -Wall
mode=1
ifeq ($(mode), 0) # production mode 
CFLAGS+=-O3 # max optimization
endif

ifeq ($(mode), 1) # development mode 
profile=1 # development mode, otherwise '0'
CFLAGS+=-O0 -g -ggdb
endif

CFLAGS+=-MMD  # generate dependency .d files
CFLAGS+=-I../include # custom 
LDLIBS=-lm # libm - basic maths
job=solus
input=strainer
srcs=$(job).c $(shell ls ../ext/*.c)
objs=$(srcs:.c=.o)
deps=$(objs:.o=.d)
ifeq ($(profile), 1)
CFLAGS+=-pg
LDFLAGS+=-pg
endif

$(job) : $(objs)
-include $(deps)
# %: %.o # all targerts for all objects
# obj/%.o: %.c # to be tested
trgs=$(job) # single target in our case 

### command line parameters:
N=3 # communicator size
C=15 # number of active centers in the data; see strainer.pl output
X=16 # eXclusion limit; currently the number of points closest to the center to retain
W=0 # number of points to weight

.PHONY: run
run :
	mpirun -n $N ./$(job) -c $C -x $X -w $W \
		-i ./$(input)-p32-FE-MN.ssv \
		-o ./$(job).ssv \
		2> ./$(job).err

.PHONY: clean
clean:
	$(RM) $(objs)
	$(RM) $(deps)

### everything below is relevant only in dev mode

# cscope.out
cscope.out: $(srcs)
	cscope -f $@ -I$(CURDIR)/inc -bq $^

# ctags
tags: $(srcs)
	ctags -f $@ -R $(CURDIR) $^

ifeq ($(profile), 1)

gmon.out: $(job)

gprof.out: gmon.out
	gprof $(CURDIR)/$(job) > $@

.PHONY: prof
prof: gprof.out
endif

.PHONY: more
more: cscope.out tags prof

.PHONY: cleanall
cleanall:
	$(RM) $(trgs)
	$(RM) $(objs)
	$(RM) $(deps)
	$(RM) $(job).ssv $(job).err
	$(RM) cscope.out cscope.out.in cscope.out.po
	$(RM) tags
ifeq ($(profile), 1)
	$(RM) gmon.out
	$(RM) gprof.out
endif
