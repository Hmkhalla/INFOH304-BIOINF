#MKHALLATI Hassan 000459680
test_db : Database.o 
	g++ Database.o test_db.cpp -o test_db
	rm -rf *.o
Database.o : Database.cpp
	g++ -c Database.cpp -o Database.o
Protein.o : Protein.cpp
	g++ -c Protein.cpp -o Protein.o
clean : 
	rm -rf *.o
