gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c hgtInit.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c ThetaOfT.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c GramAtN.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c GramNearT.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c RSbuildcoeff.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c RSremainder.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c RSmainTerm.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c HardyZcalc.c 
ar rcs libhgt.a hgtInit.o ThetaOfT.o GramAtN.o GramNearT.o RSbuildcoeff.o RSremainder.o RSmainTerm.o HardyZcalc.o

