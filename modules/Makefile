include ../make.sys

objects = constants.o para.o input.o dosdata.o wanndata.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a
