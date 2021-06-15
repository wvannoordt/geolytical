EXE := arc
##MAIN := arcLength_airfoil.cpp
MAIN := a_airfoil.cc
test: main
	./arc
main:
	g++ -I. -I${GEOLYTICAL}/include ${MAIN} -o ${EXE} -L${GEOLYTICAL}/lib -lgeolytical

clean:
	rm -f naca arc *~ *vtk
