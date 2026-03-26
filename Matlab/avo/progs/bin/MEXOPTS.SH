# mexopts.sh
	CC='g++'
	CFLAGS='-DUnix \
		-I/home/bioel/wilberth/c -I/home/bioel/wilberth/bem/c \
		-I/home/bioel/wilberth/include \
		-I/home/bioel/thom/c -I/home/bioel/thom/bem/c \
		-DLinux -Dsystemid=\"Linux\" -DMEX'
	CLIBS='-lm -lc\
		-L/home/bioel/wilberth/lib -lnm\
		-L/home/bioel/thom/lib -lThomlinux -lBemlinux'
	COPTIMFLAGS='-O -DNDEBUG'
	CDEBUGFLAGS='-g'
#
	FC='fort77'
	FFLAGS=''
	FLIBS='-lf2c -Wl,--defsym,MAIN__=mexfunction_'
	FOPTIMFLAGS='-O'
	FDEBUGFLAGS='-g'
#
	LD=$CC
	LDFLAGS='-shared'
	LDOPTIMFLAGS=''
	LDDEBUGFLAGS=''
