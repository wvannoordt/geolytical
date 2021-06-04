EXE := naca
MAIN := $(wildcard *.cc)
run: main
	./naca
main:
	g++ -I. -I${GEOLYTICAL}/include ${MAIN} -o ${EXE} -L${GEOLYTICAL}/lib -lgeolytical

clean:
	rm -f test.x *~ *vtk
