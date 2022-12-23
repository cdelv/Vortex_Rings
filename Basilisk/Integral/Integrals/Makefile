

all:
	g++ -fpic -shared Integral.cpp -o libintegral.so
	g++ -fpic -shared connector.cpp -L. -lintegral -o libconnector.so