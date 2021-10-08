OBJECTS = main.o met_microbe_coupling.o node_node_ABC.o printWs.o matrix.o

getWmatrix: $(OBJECTS)
	g++ -o getWmatrix $(OBJECTS)

main.o: main.cpp functions.h matrix.h
	g++ -c -std=c++11 main.cpp

met_microbe_coupling.o: met_microbe_coupling.cpp functions.h
	g++ -c -std=c++11 met_microbe_coupling.cpp

node_node_ABC.o: node_node_ABC.cpp functions.h matrix.h
	g++ -c -std=c++11 node_node_ABC.cpp

printWs.o: printWs.cpp functions.h matrix.h
	g++ -c -std=c++11 printWs.cpp

matrix.o: matrix.cpp matrix.h
	g++ -c -std=c++11 matrix.cpp

clean:
	-rm -f $(OBJECTS)
