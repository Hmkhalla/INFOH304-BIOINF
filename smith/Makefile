#MKHALLATI Hassan 000459680
fast=-Ofast
test_db : Database.o 
	g++ Database.o test_db.cpp -o test_db $(fast)
	rm -rf *.o
Database.o : Database.cpp
	g++ -c Database.cpp -o Database.o $(fast)
Protein.o : Protein.cpp
	g++ -c Protein.cpp -o Protein.o $(fast)
clean : 
	rm -rf *.o
