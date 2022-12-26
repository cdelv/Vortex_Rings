

all:
	g++ -O3 -fpic -shared Integral.cpp -o libintegral.so
	g++ -O3 -fpic -shared connector.cpp -L. -lintegral -o libconnector.so