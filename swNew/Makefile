fast=-Ofast
projet: Database.o Sequence.o Algorithm.o
	g++ Database.o Sequence.o Algorithm.o projet.cpp -o projet $(fast) -lpthread
Database.o: Database.cpp 
	g++ -c Database.cpp -o Database.o $(fast)
Sequence.o: Sequence.cpp
	g++ -c Sequence.cpp -o Sequence.o $(fast)
Algorithm.o: Algorithm.cpp
	g++ -c Algorithm.cpp -o Algorithm.o $(fast) -lpthread
clean: 
	rm -rf *.o
	
	
