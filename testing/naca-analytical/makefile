EXE := naca
MAIN := $(wildcard *.cc)

.PHONY: test
test: run

run: main
	./naca
main:
	g++ -I. -I${GEOLYTICAL}/include ${MAIN} -o ${EXE} -L${GEOLYTICAL}/lib -lgeolytical

clean:
	rm -f naca *~ *vtk
