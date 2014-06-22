INCLUDEFLAGS := -I$(HOME)/local/include/
CFLAGS := -g -std=c99 -Wall
CPPFLAGS := -g -Wall -std=gnu++11 -fopenmp -msse3 -O3
LFLAGS := -L$(HOME)/local/lib/ -L/usr/lib/x86_64-linux-gnu/ -g -O3
# -ltcmalloc -lprofiler
LIBS := -lgsl -lgslcblas -fopenmp -ljansson

