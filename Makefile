FC = mpif90 -O2 ## gfotran
#FC = mpiifort -O3 -xHOST -ipo -ip ## intel
#FC = mpifrtpx -O3 -Kfast ##FX100@Nagoya

LN = -llapack -lblas

VPATH = src:object
SRC = $(shell cd src ;ls *.f90 ;cd ..)
OBJ = $(SRC:.f90=.o)
OBJ_dir = $(addprefix object/,$(OBJ))

PROG = hchain

$(PROG):global_variables.o $(OBJ)
	$(FC) -o $(PROG) $(OBJ_dir) $(LN)

main.o:main.f90
	$(FC) -c $< $(LN);mv $@  object 
%.o:%.f90
	$(FC) -c $< $(LN);mv $@  object 


clean:
	rm  -f  object/*.o  *.mod ${PROG}
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod ${PROG} */#* *.out *.log
