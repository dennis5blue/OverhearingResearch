all:main.cpp
	g++ -L/usr/local/lib main.cpp -o main
clean:
	rm -f main
