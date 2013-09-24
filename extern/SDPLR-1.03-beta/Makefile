all: gsl sdplr

gsl:
	cd gsl-1.5/poly; make

sdplr:
	cd source; make

clean:
	cd gsl-1.5/poly; make clean
	cd source; make clean

cleanall: clean
	rm -rf sdplr
	rm -rf sdplr.exe
	rm -rf mexsdplr.*
	rm -rf lib/libgsl.a

mingw: mingw_gsl mingw_sdplr

mingw_gsl:
	cd gsl-1.5\poly & mingw32-make

mingw_sdplr:
	cd source & mingw32-make

mingw_clean:
	cd gsl-1.5\poly & mingw32-make mingw_clean
	cd source & mingw32-make mingw_clean

mingw_cleanall: mingw_clean
	del sdplr
	del sdplr.exe
	del mexsdplr.*
	del lib\libgsl.a
