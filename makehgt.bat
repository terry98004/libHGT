gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c hgtInit.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c ThetaOfT.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c GramAtN.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c GramNearT.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c RSbuildcoeff.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c RSremainder.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c RSmainTerm.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c HardyZcalc.c 
ar rcs libhgt.a hgtInit.o ThetaOfT.o GramAtN.o GramNearT.o RSbuildcoeff.o RSremainder.o RSmainTerm.o HardyZcalc.o

