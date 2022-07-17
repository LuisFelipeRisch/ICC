    CC     = gcc -g
    CFLAGS = 
    LFLAGS = -lm

      PROG = main
      OBJS = utils.o \
             edoHandler.o \
						 coefficients.o \
             solve.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h utils.h edoHandler.h coefficients.h solve.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)

