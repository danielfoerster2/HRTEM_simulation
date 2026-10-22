FC = mpifort
FFLAGS ?= -O2 -Wall -Wextra
#FFLAGS ?= -Wall -Wextra -g -O0 -fcheck=all -fbacktrace
LDLIBS ?= /lib/x86_64-linux-gnu/libfftw3.so.3 /lib/x86_64-linux-gnu/libopenblas.so.0

PREP_OBJS = main.o utils_io.o constants.o random_utils.o nano_process.o slc.o fft.o msa.o wavimg.o

all: main

main: $(PREP_OBJS)
	$(FC) $(FFLAGS) $^ -o $@ $(LDLIBS)

constants.o: constants.f90
	$(FC) $(FFLAGS) -c $< -o $@

utils_io.o: utils_io.f90 constants.o random_utils.o
	$(FC) $(FFLAGS) -c $< -o $@

nano_process.o: nano_process.f90 constants.o
	$(FC) $(FFLAGS) -c $< -o $@

random_utils.o: random_utils.f90 constants.o
	$(FC) $(FFLAGS) -c $< -o $@

slc.o: slc.f90 constants.o fft.o
	$(FC) $(FFLAGS) -c $< -o $@

msa.o: msa.f90 constants.o fft.o
	$(FC) $(FFLAGS) -c $< -o $@

wavimg.o: wavimg.f90 constants.o fft.o
	$(FC) $(FFLAGS) -c $< -o $@

fft.o: fft.f90 constants.o
	$(FC) $(FFLAGS) -c $< -o $@


main.o: main.f90 utils_io.o constants.o random_utils.o nano_process.o slc.o fft.o msa.o wavimg.o
	$(FC) $(FFLAGS) -c $< -o $@

run: main
	mpirun -np 1 ./main

clean:
	rm -f *.o *.mod tmp.dat particle_descriptors.dat \
	      descriptors.dat xyz_file_list.tmp records_rank_*.tmp \
		  main
		
.PHONY: all run clean
