all:main.cpp
	g++ -L/opt/boost_1_55_0/stage/lib main.cpp -o main
clean:
	rm -f main
