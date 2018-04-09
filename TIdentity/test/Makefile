SOURCE=testIdentity.C
MYPROGRAM=testIden
RINCLUDES= $(shell root-config --incdir)
RLIBRARIES= $(shell root-config --libs  | sed 's/-pthread / /')
#TIdentityLIBS=../lib/libTIdentity.so

CC = g++ -m64

all: $(MYPROGRAM)

$(MYPROGRAM): $(SOURCE)

	$(CC) $(INCLUDE) -L ../lib  -lTIdentity -I../  -I$(RINCLUDES) $(SOURCE) -o$(MYPROGRAM) $(RLIBRARIES) $(TIdentityLIBS)

clean:

	rm -f $(MYPROGRAM)
