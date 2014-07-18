INCLUDEFLAGS := -I/usr/local/include/
CFLAGS := -g -std=c99 -Wall
CPPFLAGS := -g -Wall -std=gnu++11 -msse3 -O3
LFLAGS := -L/usr/local/lib/ -g -O3
# -ltcmalloc -lprofiler
LIBS := -lgsl -lgslcblas -ljansson
