all:
	g++ -std=c++23 bitpal.cpp -O3 -o bitpal

clean:
	rm -f bitpal