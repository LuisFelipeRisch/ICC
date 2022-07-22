    CC     = gcc -g

    CFLAGS = -I/usr/include
    LFLAGS = -lm -L/usr/lib -lmatheval

      PROG = labZero
      OBJS = utils.o \
             ZeroFuncao.o \
             functions.o

.PHONY: clear purge clean purge distclean all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) $(PROG).o
	$(CC) -o $@ $^ $(LFLAGS)

clear:
	@rm -f *~ *.bak

purge:   clear
	@rm -f *.o core a.out
	@rm -f $(PROG)
