include ../make.sys

objects = flib.o interpolate_bands.o

LIBS = ../modules/mod.a

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

wanndos: wanndos.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o wanndos.x wanndos.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

clean:
	rm -rf *.o *.mod *.x
