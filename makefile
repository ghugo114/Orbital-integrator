CC = gcc
CFLAGS = -O2


LIBS = /lib/x86_64-linux-gnu/libm.so.6 
LIBDIRS = /xgrafix_instal/lib/libXGC.a
LIBDIRS = /usr/local/lib/
LIBDIRS = /xgrafix
# Link the two object files into an executable file "planets".



select_trayect: select_trayect.o  
	cc select_trayect.o -o select_trayect $(LIBS)


# Compile source file "planets.c" to make an object file "planets.o".

select_trayect.o: select_trayect.c planets.h 
	cc -c select_trayect.c 




