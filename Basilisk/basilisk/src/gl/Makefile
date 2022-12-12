-include ../config

CFLAGS += -O2 -I..

all: libglutils.a libfb_dumb.a libfb_osmesa.a libfb_glx.a

libglutils.a: trackball.o utils.o polygonize.o \
		og_stroke_mono_roman.o parser.o
	ar cr $@ $^

libfb_dumb.a: fb_dumb.o
	ar cr $@ $^

fb_dumb.o: gl.h

libfb_osmesa.a: fb_osmesa.o
	ar cr $@ $^

libfb_glx.a: fb_glx.o OffscreenContextGLX.o fbo.o
	ar cr $@ $^

utils.o: utils.h

parser.tab.c: parser.y
	bison parser.y

parser: parser.tab.c parser.h
	$(CC) $(CFLAGS) -DSTANDALONE=1 -o parser parser.tab.c -lm

parser.o: parser.tab.c parser.h
	$(CC) $(CFLAGS) -c -o parser.o parser.tab.c

clean:
	rm -f *.o *.a
