include ../Makefile.inc

LIB_DIRS_2 = -L../lib

APPL_SRCS = \
          copystructures.c \
          dataoper.c \
          eigen.c \
          initialize.c \
          lbfgs.c \
          linesearch.c \
          main.c \
          misc.c \
          params.c \
          rankreduce.c \
          readdata.c \
          sdplrlib.c \
          timefuncs.c \
          util.c

APPL_OBJS = $(APPL_SRCS:.c=.o)

sdplr: $(APPL_OBJS) 
	$(CC) $(CFLAGS) $(LIB_DIRS) $(LIB_DIRS_2) -o ../sdplr $(APPL_OBJS) $(LIBS)

$(APPL_OBJS): $(OBJECT)%.o : $(SOURCE)%.c
	$(CC) $(CFLAGS) $(DEFINES) -c $< -o $@

clean:
	rm -f *.o

mingw_clean:
	del *.o
