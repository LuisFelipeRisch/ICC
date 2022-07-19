    CC     = gcc -g
    CFLAGS = 
    LFLAGS = -lm

      PROG = main
      OBJS = utils.o \
             solve.o \
             edoHandler.o \
						 coefficients.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h utils.h solve.h edoHandler.h coefficients.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)

