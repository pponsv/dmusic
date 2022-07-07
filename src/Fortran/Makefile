build:
	gfortran -Wno-tabs -c helper.f90
	f2py --f90flags="-Wno-tabs -O2" -llapack -c main.f90 helper.f90 -m dmusic_f

run:
	python3 test.py

all: build run

clean:
	rm -r *.cpython* *.o *.mod

build_debug:
	gfortran -Wno-tabs -c helper.f90
	f2py --f90flags="-Wno-tabs -O2" -llapack -c debug.f90 helper.f90 -m debug

run_main:
	python3 test_main.py

all_main: build_main run_main
# main:
# 	gfortran -Wno-tabs -c helper.f90
# 	gfortran -o main.f90 helper.f90 -m main
