ssnew.x: besselnew.o square_routines.o squarenew.o GOE.o rgnf_lux.o matrix_stuff.o ~/bin/minpack.o
	gfortran GOE.o rgnf_lux.o matrix_stuff.o besselnew.o square_routines.o ~/bin/minpack.o squarenew.o -fdefault-real-8 -fdefault-double-8 -o ssnew.x -L /Users/mehtan/Code/ARPACK/ARPACK/ -larpack_OSX -L/usr/local/opt/lapack/lib/ -llapack -lblas

squarenew.o: squarenew.f90
	gfortran -fcheck-new -fcheck=bounds -Wargument-mismatch -Winteger-division -Wsurprising -Wintrinsic-shadow -fdefault-real-8 -fdefault-double-8 -c squarenew.f90

besselnew.o: besselnew.f
	gfortran -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -fno-range-check -c besselnew.f

square_routines.o: square_routines.f
	gfortran -fdefault-real-8 -fdefault-double-8 -c square_routines.f

GOE.o:	GOE.f
	gfortran  -ffixed-line-length-132 -c GOE.f

rgnf_lux.o: rgnf_lux.f
	gfortran  -ffixed-line-length-132 -c rgnf_lux.f

matrix_stuff.o:	matrix_stuff.f
	gfortran  -ffixed-line-length-132 -c matrix_stuff.f

~/bin/minpack.o: ~/bin/minpack.f90
	gfortran    -fdefault-real-8 -fdefault-double-8 -c ~/bin/minpack.f90
