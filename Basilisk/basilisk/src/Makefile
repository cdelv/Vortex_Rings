CFLAGS += -O2

# these are not Basilisk programs
EXCLUDE = qcc.c include.c postproc.c rotate.c stencils.c qplot.c

all: libast qcc qplot libkdt literatec libfb_dumb libws bview2D bview3D
	@chmod +x ppm2mpeg ppm2mp4 ppm2ogv ppm2gif runtest \
		sequence page2html
	@test -f xyz2kdt || ln -s kdt/xyz2kdt
	@test -f kdtquery || ln -s kdt/kdtquery

libast:
	cd ast && make

libkdt:
	cd kdt && make

literatec:
	cd darcsit && make && cd cgi-bin && make

qcc: qcc.c include.o postproc.o ast/libast.a config
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" \
		-DCC99="\"$(CC99)\"" \
		-DCPP99="\"$(CPP99)\"" \
		-DCADNACC="\"$(CADNACC)\"" \
		-DBASILISK="\"$(BASILISK)\"" \
		qcc.c include.o postproc.o -o qcc -Last -last

include.o: include.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC) $(CFLAGS) -c qplot.c

qplot: qplot.o
	$(CC) $(CFLAGS) qplot.o -o qplot

draw_get.h: draw.h params.awk
	awk -f params.awk < draw.h > draw_get.h

draw_json.h: draw.h json.awk
	awk -f json.awk < draw.h > draw_json.h

include.c: include.lex
	flex -P inc -o include.c include.lex

libfb_dumb:
	cd gl && make libfb_dumb.a libglutils.a

libws:
	cd wsServer && make

postproc.c: postproc.lex
	flex -P post -o postproc.c postproc.lex

bview2D: qcc bview.c draw_get.h draw_json.h bview.s
	qcc $(CFLAGS) -DDUMBGL -autolink bview.c -o bview2D -lfb_dumb -lm

bview3D: qcc bview.c draw_get.h draw_json.h bview.s
	qcc $(CFLAGS) -DDUMBGL -autolink -grid=octree bview.c -o bview3D -lfb_dumb -lm

alltags:
	cd navier-stokes && make tags
	cd layered && make tags
	cd ehd && make tags
	cd examples && make tags
	cd test && make tags
	make tags
	cd navier-stokes && make itags
	cd layered && make itags
	cd ehd && make itags
	cd examples && make itags
	cd test && make itags
	make itags

etags:
	etags *.h grid/*.h

checklinks:
	$(LINKCHECKER) 	$(BASILISK_URL)/src/README 		\
			$(BASILISK_URL)/src/test/README 	\
			$(BASILISK_URL)/src/examples/README | 	\
		tee checklinks.log

checklinksfast:
	wget --spider -nd -nv -r \
		--reject-regex '.*[?]changes=.*' \
		--reject-regex '.*[?]history' \
		$(BASILISK_URL) 2>&1 | \
		grep -v ^unlink: | tee checklinks.log

changelog:
	darcs changes > ChangeLog

dist:
	darcs dist

diff:
	cd .. && tar czvf src/diff.tgz `darcs whatsnew -s | \
		sed 's/. .\/.*\/$$//g' | awk '{print $$2}'`

include $(BASILISK)/Makefile.defs
