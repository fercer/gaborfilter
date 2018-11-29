SRCDIR:=./include
LIBDIR:=./lib
ODIR:=./obj

EXAMPLESDIR:=./examples
EXAMPLESSRCDIR:=$(EXAMPLESDIR)/src
EXAMPLESODIR:=$(EXAMPLESDIR)/obj
EXAMPLESBINDIR:=$(EXAMPLESDIR)/bin

SRCS:=$(wildcard $(SRCDIR)/*.c)
HDRS:=$(wildcard $(SRCDIR)/*.h)
OBJS:=$(patsubst $(SRCDIR)/%,$(ODIR)/%, $(patsubst %.c,%.o, $(SRCS)))

LIBS:=libgabor.so

make: $(OBJS) $(LIBDIR)
	gcc -shared -L/usr/lib/x86_64-linux-gnu -lfftw3 -o $(LIBDIR)/$(LIBS)

$(OBJS): $(ODIR)
	gcc -c -fpic -I/usr/include/ $(SRCS) -o $(OBJS) -std=c99

$(ODIR):
	mkdir $(ODIR)

$(LIBDIR):
	mkdir $(LIBDIR)

examples: $(EXAMPLESOBJS)
	cp $(LIBDIR)/$(LIBS) $(EXAMPLESDIR)/$(LIBS)

example_gabor: $(EXAMPLESBDIR) $(EXAMPLESODIR)/test_gabor.bin
	gcc $(EXAMPLESODIR)/test_gabor.bin -c $(EXAMPLESSRC)/test_gabor.c
